"""Explicit systematic component encoders, table values, and RAM program."""
from hashlib import sha256

LAYOUT = (30,) * 33 + (10,)
MAIN_COLUMNS = tuple(
    int.from_bytes(sha256(f"wide-product-a:30:{i}".encode("ascii")).digest()[:4], "little")
    & ((1 << 30) - 1)
    for i in range(30)
)
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
    result = 0
    offset = 0
    for k in LAYOUT:
        result |= pack_block((message >> offset) & cap(k), k) << (2 * offset)
        offset += k
    return result


def table_addresses():
    return {30: (1 << 31, 1 << 30), 10: (1 << 11, 1 << 10)}, 3 << 30


def schedule():
    for j in range(32):
        yield "load", f"Y{j}", j
    for j in range(16):
        yield "move", f"U{j}", 0
    addresses, _ = table_addresses()
    offset = 0
    for k in LAYOUT:
        index, shift = divmod(2 * offset, 64)
        if shift:
            yield "shr", "R", f"Y{index}", shift
        else:
            yield "move", "R", f"Y{index}"
        if shift + 2 * k > 64:
            yield "shl", "T", f"Y{index + 1}", 64 - shift
            yield "or_", "R", "R", "T"
        yield "and_", "A", "R", cap(k)
        yield "shr", "P", "R", k
        residue_bits = (64 - shift - k) if shift + 2 * k <= 64 else (64 - k)
        if residue_bits != k:
            yield "and_", "P", "P", cap(k)
        p_base, _ = addresses[k]
        yield "add", "T", "A", p_base
        yield "load", "T", "T"
        yield "xor", "P", "P", "T"
        yield "load", "T", "P"
        yield "xor", "A", "A", "T"
        index, shift = divmod(offset, 64)
        if shift:
            yield "shl", "T", "A", shift
            yield "or_", f"U{index}", f"U{index}", "T"
        else:
            yield "move", f"U{index}", "A"
        if shift + k > 64:
            yield "shr", f"U{index + 1}", "A", 64 - shift
        offset += k
    for j in range(16):
        yield "store", None, 32 + j, f"U{j}"
