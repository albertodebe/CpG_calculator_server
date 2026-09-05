# my_server.py
import socket
import threading
import json
import argparse
import requests
from socket import AF_INET, SOCK_STREAM


# Gets input file in one comapct line  and remove header, while removing fasta header. 
def preprocess_file(data):
    with open(data) as f:
        file = f.read()
        #check if file type is fasta and if so remove the header
        if file[0] == ">":
            lines = file.strip().splitlines()
            #make sure everything is in uppercase
            sequence = ''.join(lines[1:]).upper()
            for i in sequence:
                if i not in "ATGC":
                    raise ValueError("INVALID INPUT: sequence contains characters that are not bases")
            return sequence
        else:
            lines = file.strip().splitlines()
            sequence = ''.join(lines).upper()
            for i in sequence:
                if i not in "ATGC":
                    raise ValueError("INVALID INPUT: sequence contains characters that are not bases")
            return sequence

# Returns percentage of CpG
def calculate_cpg_percentage(seq):
    counter=0
    for i in range(len(seq)-1):
        if seq[i]=="C" and seq[i+1]=="G":
            #the counter goes up by two so every time a C is followed by a G, both are counted, otherwise in
            #a hypotetical sequence of only CGCG, the percentage would only be 50
            counter+=2
    return (counter*100)/(len(seq)-1) 


# Divides a DNA sequence into overlapping windows and calculates the CpG percentage for each window. 
# Returns a list of tuples containing the window positions and corresponding CpG content
def sliding_window_cpg(seq, window_size,step_size):
    results=[]
    for i in range(0,len(seq)-window_size+1,step_size):
        window=(i, i+window_size)
        window_seq=seq[i:i + window_size]
        percent= calculate_cpg_percentage(window_seq)
        results.append((window[0],window[1], percent))
    return results


#returns windows with over 60% CpG 
def detect_cpg_islands(results):
    return [elem for elem in results if elem[2] >= 60.0]


## with cCREs I can't directly query the database, I have to download the BED files, so insted I'm using the ENSEMBLE
## REST API to directly retrieve the the annotation. I assume that the sequence has already been
## aligned (possibly through BLAST) and so the genomic cordinates are already known

#Finds the genetic features overlapping with the given cordinates
def query_ensembl_gene_region(chrom, start, end):
    url = f"https://rest.ensembl.org/overlap/region/human/{chrom}:{start}-{end}?feature=gene;feature=regulatory"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)
    if response.ok:
        #returns a list of features that overlap with our sequence:
        #[
        # { 
        #  "id": "ENSG00000186092",
        #  "feature_type": "gene",
        #  "external_name": "OR4F5",
        #  "start": 11869,
        #  "end": 14412, 
        #   ...
        # }, 
        # ...
        #]
        return response.json()
    else:
        return [{"annotation": "API error or no gene found"}]


def handle_client(conn):
    try:
        #sets an upper bound for the request, also decodes (converts intro a string that can be processed as text)
        request = conn.recv(65536).decode() 
        data = json.loads(request)
        sequence = preprocess_file(data['sequence'])
        response = {}
        # differentiate between global CpG count and sliding Window CpG profile 
        if data['mode'] == 'global':
            response['cpg_percentage'] = calculate_cpg_percentage(sequence)
        elif data['mode'] == 'sliding':
            win = data['window_size']
            step = data['step_size']
            chrom=data['chrom']
            start=data['start']
            end=data['end']
            profile = sliding_window_cpg(sequence, win, step)
            cpgs = detect_cpg_islands(profile)
            #in case there are no gpg islands and annotation is empty
            if cpgs:
                annotations = query_ensembl_gene_region(chrom, start, end)
            else:
                annotations = [{"annotation": "No CpG islands detected; annotation skipped"}]
            response['profile'] = profile
            # summary statistics of the CpG content in the sliding windows
            response['summary_stats'] = {
                "max": max(p[2] for p in profile),
                "min": min(p[2] for p in profile),
                "mean": sum(p[2] for p in profile) / len(profile)
            }
            response['cpg_islands'] = cpgs
            response['annotations'] = annotations
        # converts the response dictionary into a JSON and ecodes it    
        conn.sendall(json.dumps(response).encode())
    except Exception as e:
        conn.sendall(json.dumps({"error": str(e)}).encode())
    finally:
        conn.close()


# Starts a server that listens on a specified host and port and accepts incoming clients creating separate threads
def start_server(host, port):
    s = socket.socket(family=AF_INET, type=SOCK_STREAM, proto=0)
    s.bind((host, port))
    s.listen(5)
    print(f"[Server] Listening on {host}:{port}")
    while True:
        conn, addr = s.accept()
        # creates a new thread for every client and executes handle_client
        thread = threading.Thread(target=handle_client, args=(conn,))
        thread.daemon = True
        thread.start()


if __name__ == "__main__":
    #creates an argument parser object, which allows to pass arguments from the command line
    parser = argparse.ArgumentParser()
    #127.0.0.1 allows the computer to comunicate with itself
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=9999)
    args = parser.parse_args()
    start_server(args.host, args.port)