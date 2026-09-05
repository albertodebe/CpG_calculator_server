# CpG_calculator_server

A multithreaded TCP socket server implemented in Python that focuces on CpG content evaluation.
After preprocessing it calculates the CpG percentage, either in global mode (ie for the whole sequence) or for overlapping windows, whose length and step size can be personalized. GpG islands (above 60%) are then queried against the ENSEMBL API.
