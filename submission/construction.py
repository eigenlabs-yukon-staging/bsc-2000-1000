"""Explicit systematic component encoders, table values, and RAM program."""
LAYOUT = (30,) * 33 + (10,)
GENERATOR = 14978968423


def polynomial_remainder(value):
    while value.bit_length() > 33:
        value ^= GENERATOR << (value.bit_length() - 34)
    return value


PUNCTURED = (0,7,26)
KEPT = tuple(j for j in range(33) if j not in PUNCTURED)


def main_column(i):
    remainder = polynomial_remainder(1 << (33 + i))
    return sum(((remainder >> j) & 1) << b for b, j in enumerate(KEPT))


MAIN_COLUMNS = tuple(main_column(i) for i in range(30))
TAIL_COLUMNS = (980, 731, 346, 811, 541, 494, 407, 357, 441, 744)


def parity(data, columns):
    result = 0
    for i, column in enumerate(columns):
        if (data >> i) & 1:
            result ^= column
    return result


def cap(k):
    return (1 << k) - 1


def block_columns(k):
    return MAIN_COLUMNS if k == 30 else TAIL_COLUMNS


def syndrome_columns(k):
    return block_columns(k) + tuple(1 << i for i in range(k))


def pack_block(data, k):
    return data | (parity(data, block_columns(k)) << k)


def encode(message):
    result=message
    offset=0
    for k in LAYOUT:
        data=(message>>offset)&cap(k)
        result |= parity(data,block_columns(k))<<(1000+offset)
        offset+=k
    return result


def table_addresses():
    return {30: (1 << 31, 1 << 30), 10: (1 << 11, 1 << 10)}, 3 << 30


def extract(destination,start,width):
    index,shift=divmod(start,64)
    if shift:
        yield "shr",destination,f"Y{index}",shift
    else:
        yield "move",destination,f"Y{index}"
    if shift+width>64:
        yield "shl","T",f"Y{index+1}",64-shift
        yield "or_",destination,destination,"T"
    if shift+width!=64:
        yield "and_",destination,destination,cap(width)


def schedule():
    for j in range(32):
        yield "load",f"Y{j}",j
    for j in range(16):
        yield "move",f"U{j}",f"Y{j}"
    addresses,_=table_addresses()
    offset=0
    for k in LAYOUT:
        yield from extract("A",offset,k)
        yield from extract("P",1000+offset,k)
        p_base,_=addresses[k]
        yield "add","T","A",p_base
        yield "load","T","T"
        yield "xor","P","P","T"
        yield "load","R","P"
        index,shift=divmod(offset,64)
        output=f"U{index}"
        if shift:
            yield "shl","T","R",shift
            yield "xor",output,output,"T"
        else:
            yield "xor",output,output,"R"
        if shift+k>64:
            next_output=f"U{index+1}"
            yield "shr","T","R",64-shift
            yield "xor",next_output,next_output,"T"
        offset+=k
    yield "and_","U15","U15",cap(40)
    for j in range(16):
        yield "store",None,32+j,f"U{j}"
