import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input', help = 'input file')
args = parser.parse_args()

def checkDuplicates(listOfElems):
    '''
    Checks if there are duplicate elements in a given list.
    '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

def main():
    '''
    Reads the input gene order file, checks for duplicates, 
    and print out number of total entries and number of unique entries if duplicates exist.
    '''
    with open(args.input,'r') as f:
        for line in f:
            orf_number_list = line.split(' ')
    if checkDuplicates(orf_number_list) == False:
        print('No duplicate entries')
        print(len(orf_number_list))
    else:
        print('Contains duplicate entries')
        print(len(orf_number_list))
        print(len(set(orf_number_list)))

if __name__ == '__main__':
    main()
