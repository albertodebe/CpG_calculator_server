import socket
import json
import argparse


# Sends a JSON request and returns server's response, which gets decoded
def send_request_to_server(host, port, input_data):
    s=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((host, port))
    #converts the input data into into a JSON string
    s.sendall(json.dumps(input_data).encode())
    #max file size is set to 10MB
    response = s.recv(10_485_760).decode()
    return json.loads(response)


# Creates a parser object that can take command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1",help="server IP address")
    parser.add_argument("--port", type=int, default=9999,help="port number")
    parser.add_argument("--file", required=True, help="path/to/fasta(or the raw sequence)")
    parser.add_argument("--mode", required=True, choices=["global", "sliding"],help="anaylis mode: Global for the whole sequence, sliding for sliding door method")
    parser.add_argument("--window", type=int,help="window size for sliding window")
    parser.add_argument("--step", type=int,help="step size for sliding window")
    parser.add_argument("--chrom",help="chromosome number")
    parser.add_argument("--start", type=int,help="start coordinate on the chromosome")
    parser.add_argument("--end", type=int,help="end coordinate on the chromosome")
    args = parser.parse_args()

    if args.mode=="global":
        input_data = {"sequence": args.file,
                      "mode": args.mode,
                    }
    elif args.mode=="sliding": 
        input_data = {"sequence": args.file,
                      "mode": args.mode,
                      "window_size": args.window,
                      "step_size": args.step,
                      "chrom": args.chrom,
                      "start": args.start,
                      "end": args.end
                     }
    else:
        raise Exception("mode entered is not valid, select either 'global' or 'sliding' ")


    response = send_request_to_server(args.host, args.port, input_data)
    # Converts the python dictionary into a JSON string
    print(json.dumps(response, indent=2))
