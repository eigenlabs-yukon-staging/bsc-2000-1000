"""Fixed row identities, adaptive cap tables, and their exact certificates."""
from bisect import bisect_right
from collections import Counter
from functools import lru_cache
import json
from math import comb
from pathlib import Path

from construction import (BUDGETS, BYTES, CAP_OFFSETS, GROUPS, LIMIT, M,
    NIBBLES, PARITY, PROFILE_BASE, PROFILE_COUNTS, WIDTHS, WORDS, cost)
from qr import N, ORBITS, ROW_CAP, error_pattern, orbit_labels

SCALE = 10**60
SCHEDULE = WIDTHS*(GROUPS-1)


@lru_cache(None)
def small_values(width):
    weights = tuple(i.bit_count() for i in range(1<<width))
    parity = tuple(PARITY[width](i) for i in range(1<<width))
    return weights, parity


@lru_cache(None)
def small_profile(width, syndrome):
    weights, parity = small_values(width)
    sizes = [0]*(2*width+1)
    for weight, image in zip(weights,parity,strict=True):
        sizes[weight+weights[image^syndrome]] += 1
    return tuple(sizes)


@lru_cache(None)
def data():
    stored = json.loads((Path(__file__).resolve().parent/"profiles.json").read_text())
    qr = {tuple(row["h"]): row["multiplicity"] for row in stored}
    histograms = {24:qr}
    for width in (12,4):
        histograms[width] = Counter(small_profile(width,s) for s in range(1<<width))
    for width, histogram in histograms.items():
        assert len(histogram) == PROFILE_COUNTS[width]
        assert sum(histogram.values()) == 1<<width
    return [(width, dict(sorted(histograms[width].items()))) for width in (24,12,4)]


@lru_cache(None)
def profile_index(width, index):
    if width == 24:
        representative = orbit_labels()[1][index]
        sizes = [0]*49
        for error in range(N):
            sizes[error_pattern(error,representative).bit_count()] += 1
        histogram = tuple(sizes)
    else:
        histogram = small_profile(width,index)
    types = dict(data())[width]
    return tuple(types).index(histogram)


def kept_counts(width, histogram):
    remaining = ROW_CAP if width==24 else 1<<width
    for total in histogram:
        kept = min(total,remaining)
        remaining -= kept
        yield kept


def statistics(den):
    result = {}
    for width, histogram in data():
        by_type = []
        nominal = [9**(2*width-w) for w in range(2*width+1)]
        tilt = [(den-1)**(2*width-w) for w in range(2*width+1)]
        for profile, multiplicity in histogram.items():
            normalizer = sum(a*b for a,b in zip(profile,nominal,strict=True))
            counts, masses = Counter(), Counter()
            for kept, weight, likelihood in zip(kept_counts(width,profile),nominal,tilt,strict=True):
                if kept:
                    price = cost(normalizer,weight)
                    counts[price] += kept
                    masses[price] += multiplicity*kept*likelihood
            by_type.append((tuple(sorted(counts.items())),tuple(sorted(masses.items()))))
        result[width] = by_type
    return result


@lru_cache(None)
def derive(den, store_caps=False):
    cube_shift = 3*WORDS+3
    square = BYTES*(5+3*WORDS)+WORDS
    multiply = (3*WORDS+3*(cube_shift+WORDS)+44*WORDS
                +NIBBLES*(5+cube_shift+3*WORDS)+WORDS)
    cube_ops = square+multiply
    types = statistics(den)
    bound = [0]*(LIMIT+1)
    success = [SCALE]*(LIMIT+1)
    caps = [None]*len(SCHEDULE)
    for index in range(len(SCHEDULE)-1,-1,-1):
        width = SCHEDULE[index]
        reserve = len(SCHEDULE)-1-index
        step = 80 if width==24 else 28
        if reserve==0:
            step += cube_ops+320
        elif (index+1)%len(WIDTHS)==0:
            step += cube_ops+180
        setup = (60 if width==24 else 40)+7

        continuation = [w+step for w in bound]
        usages = []
        for counts,_ in types[width]:
            usage = [0]*(LIMIT+1)
            for price, count in counts:
                for cap in range(price+reserve,LIMIT+1):
                    usage[cap] += count*continuation[cap-price]
            assert all(x<=y for x,y in zip(usage,usage[1:]))
            usages.append(usage)
        worst = list(map(max,zip(*usages)))
        bound = [setup+b for b in worst]

        numerator = [0]*(LIMIT+1)
        if store_caps:
            caps[index] = []
        for usage,(_,masses) in zip(usages,types[width],strict=True):
            effective = [bisect_right(usage,b)-1 for b in worst]
            assert all(cap>=r for r,cap in enumerate(effective))
            if store_caps:
                caps[index].append(effective)
            by_cap = [0]*(LIMIT+1)
            for price, mass in masses:
                for cap in range(price+reserve,LIMIT+1):
                    by_cap[cap] += mass*success[cap-price]
            for r,cap in enumerate(effective):
                numerator[r] += by_cap[cap]
        normalizer = den**(2*width)
        success = [n//normalizer for n in numerator]
        assert min(success)>=0 and max(success)<=SCALE
    work = 1000+sum(comb(M,d)*(44+560+cube_ops+bound[cap])
                    for d,cap in enumerate(BUDGETS))
    return work,caps,success


def fixed_tails(success):
    return [SCALE-success[cap] for cap in BUDGETS]


def adaptive_word(address):
    for width,length in ((24,ORBITS),(12,4096),(4,16)):
        offset = address-PROFILE_BASE[width]
        if 0<=offset<length:
            return profile_index(width,offset)
    index = bisect_right(CAP_OFFSETS,address)-1
    if 0<=index<len(SCHEDULE):
        profile,remaining = divmod(address-CAP_OFFSETS[index],LIMIT+1)
        return derive(10,True)[1][index][profile][remaining]
    raise ValueError("not an adaptive table address")
