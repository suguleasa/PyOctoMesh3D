#/usr/bin/python
# developed by Zachary Smith

def sign(number):
    return number/abs(number)

def bND(*args):
    """N-Dimensional Bresenham Line Algorithm

    Take m points in n-dimensions and return a list of integer-valued
    points along the line segments between them as given by Bresenham's
    algorithm.

    args = [p1, p2, ..., pm]
    pi = [xi0, xi1, ..., xin]
    """

    point_list = []

    for p in range(1, len(args)):
        p1 = args[p-1]
        p2 = args[p]
    
        if len(p1) != len(p2):
            print "Error: Points are different dimensions!"
            return "Error"

        # Delta
        D = [b-a for (a,b) in zip(p1, p2)]

        # Determining the index of the independent variable
        Dabs = [abs(b-a) for (a,b) in zip(p1, p2)]
        ind = Dabs.index(max(Dabs))
        del Dabs

        # The indices of the dependent variables
        dep = range(0, len(D))
        del dep[ind]

        line_segment = []
        append_me = [None]*len(p1)
        d = p1[:]
        d[ind] = None
        eps = [0]*len(d)
    
        # Stepping through along the independent variable
        for i in range(p1[ind], p2[ind], sign(D[ind])):

            append_me[ind] = i

            for j in dep:
                append_me[j] = d[j]
                eps[j] += D[j]

                # If 2*eps[j] exceeds D[ind] in magnitude, increment the
                # corresponding dependent variable and reset eps[j]
                if (abs(eps[j]) << 1) >= abs(D[ind]):
                    d[j] += sign(D[j])
                    eps[j] -= sign(eps[j])*sign(D[ind])*D[ind]

            line_segment.append(append_me[:])

        point_list.extend(line_segment)        

    point_list.append(args[-1])
        
    return point_list
