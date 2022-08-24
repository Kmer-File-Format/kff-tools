import sys


def filter_file(input_filename, output_filename):
    input_file = open(input_filename, "r")
    output_file = open(output_filename, "w")
    minimizer = ""
    cmpt = 0
    for line in input_file:
        if line.find("Minimizer") >= 0:
            minimizer = line[line.find("Minimizer")+len("Minimizer: "):].replace('\n', '')
        else:
            idx = line.find('|')
            if idx >= 0 and line.find('#') < 0:
                output_file.write((line[0:idx] + minimizer + line[idx + 1:] + '\t1').replace('\n', '') + '\n')
            else:
                continue
    return


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("ERROR")
        exit(1)
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    filter_file(input_path, output_path)
