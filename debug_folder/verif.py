import sys


if __name__ == "__main__":
    input_file = open(sys.argv[1], "r")
    cmpt = 0
    for line in input_file:
        if line.find("*") >= 0:
            cmpt += int(line[-3:-1])
    print(cmpt)
