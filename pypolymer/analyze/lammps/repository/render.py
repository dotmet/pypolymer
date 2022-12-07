import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Render picture')
    parser.add_argument(
        '-i',
        type=str,
        default='vmd.dat',
        help='input file name')
    parser.add_argument(
        '-o',
        type=str,
        default='out.png',
        help='output file name')

    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    soft_path = "E:\\Vmd\\tachyon_WIN32.exe"
    input_file = args.i
    _args = r'-aasamples 1000 -mediumshade -trans_vmd -res 2560 1440 -format BMP -o'
    output_file = args.o
    comand = soft_path+' '+input_file+' '+_args+' '+output_file
    os.system(comand)

if __name__ == '__main__':
    main()
