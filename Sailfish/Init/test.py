#!/usr/bin/python
import getopt
import sys

def main(argv=None):
    try:                                
        opts, args = getopt.getopt(argv, '', ['force='])
    except getopt.GetoptError:           
        usage()                          
        sys.exit(2)                     
    
    global force
    flag=False
    for opt, arg in opts:
        if opt == "--force":
            force=arg
            flag=True
    if not flag:
        force=6e-6
    
    print force

if __name__ == "__main__":
    main(sys.argv[1:])
