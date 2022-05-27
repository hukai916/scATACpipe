def closest_number(n, m):
    """
    Find the closest number to n that is no less than n and also divisible by m.
    Ref: https://www.geeksforgeeks.org/find-number-closest-n-divisible-m/
    """
    q = int(n / m)
    n1 = m * q
    if((n * m) > 0) :
        n2 = (m * (q + 1))
    else :
        n2 = (m * (q - 1))

    if (abs(n - n1) < abs(n - n2)) :
        if n1 < n:
            n1 += m
        return n1
    if n2 < n:
        n2 += m
    return n2



for i in range(1, 100):
    print(i, closest_number(i, 4))
