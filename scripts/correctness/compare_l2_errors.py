import numpy as np
import re

if __name__=='__main__':
    build_types=['serial','ompi']
    error = {}

    #floating point regular expression taken from https://stackoverflow.com/a/4703508
    error_pattern = re.compile('L2 error: [+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?')

    print 'L2 Errors for each build type:'
    for bt in build_types:
        with open(bt+'.out') as f:

            any_match=False

            for line in f:
                match = re.search(error_pattern,line)
                if match:
                    any_match=True
                    print 'L2 error for '+bt+': '+match.group(1)
                    error[bt] = np.float64(match.group(1))

            if not any_match:
                print 'ERROR!!! No L2 error found for '+bt+' build.'
                #print out cout from run with no L2 error
                with open(bt+'.out') as ff:
                    for ll in ff:
                        print ll

                exit(1)

    max_error = max([ error['serial'], error['ompi'] ])
    tol = np.finfo(np.float64).eps*max_error*100 #~10^-14
    if abs( error['serial'] - error['ompi'] ) < tol:
       exit(0)
    else:
        print 'ERROR!!! L2 Errors do not match to machine precision'
        exit(1)