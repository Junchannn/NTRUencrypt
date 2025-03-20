import numpy as np
from poly import Poly
from Crypto.Util.number import bytes_to_long, long_to_bytes
import random 

N = 107
p = 3
q = 64
df = 14
dg = 12
d = 5
    
def random_ternary_list(d1, d2, degree):
    if d1 + d2 > degree + 1:
        raise ValueError("d1 + d2 cannot exceed degree + 1 (polynomial length)")

    a = [0] * (degree + 1)

    # Choose unique indices for 1s and -1s
    ones_indices = random.sample(range(degree + 1), d1)
    remaining_indices = list(set(range(degree + 1)) - set(ones_indices))
    minus_ones_indices = random.sample(remaining_indices, d2)

    # Assign values to selected indices
    for i in ones_indices:
        a[i] = 1
    for i in minus_ones_indices:
        a[i] = -1

    return a


def random_poly_list(degree, bound):
    return np.random.randint(-bound, bound + 1, size=degree + 1).tolist()


def message_2_poly(message, n):
 
    temp = bytes_to_long(message)
    coeffs = []
    while temp > 0:
        coeffs.append(temp % n)
        temp //= n
    
    return Poly(n, coeffs)

def poly_2_message(poly):
    modulus = int(poly.n)
    temp = poly.coeffs
    message = 0

    for i in range(len(temp) - 1, -1, -1):
        message *= modulus
        message += int(temp[i])
    
    return long_to_bytes(int(message))

def complex_array_2_polys(m_list : np.ndarray[np.complex128], R : np.int64, block_size: int) -> np.ndarray[Poly]:
    
    
    poly_list = []
    for k in range(int(np.ceil(len(m_list)/ block_size))):
        temp = m_list[block_size*k:block_size*(k + 1)]
          
        coeff = []
        real_part = list(map(round,temp.real))
        imag_part = list(map(round,temp.imag))
        for i in range(len(temp)):
            coeff.append(real_part[i])
            coeff.append(imag_part[i])
            
        poly_list.append(Poly(R,coeff))
        

    return np.array(poly_list, dtype=object)

def polys_2_complex_array(poly_list: np.ndarray[np.complex128], pad = None) -> np.ndarray[np.complex128]:
    complex_array = np.array([])

    for j, p in enumerate(poly_list):
        real_part = []
        imag_part = []
        
        if ((len(p.coeffs) & 1) == 1): p.coeffs = np.concatenate((p.coeffs, [0]))            

        for i, c in enumerate(p.coeffs):
            if ((i & 1) == 0):
                real_part.append(c)
            else:
                imag_part.append(c)

        real_part = np.array(real_part, dtype=np.int64)
        imag_part = np.array(imag_part, dtype=np.int64)
        
        if pad != None:
            complex_array = np.concatenate([complex_array, pad(real_part + 1j * imag_part)])
        else: complex_array = np.concatenate([complex_array, real_part + 1j * imag_part])



    return np.array(complex_array, dtype=np.complex128)

            
        