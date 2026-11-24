"""Exact code-specific error and straight-line operation checks."""
import hashlib,json,random,sys,time
from fractions import Fraction
from pathlib import Path
HERE=next(d for d in (Path('/snapshot/files'),Path(__file__).resolve().parent) if (d/'opmodel.py').is_file())
sys.path.insert(0,str(HERE))
from construction import P48,P32,T48,T32,LUT48,LUT32,END,leaders,syndrome_luts,encode,encode_block
from decoder import PROGRAM,W_ACTUAL,decode,TableMemory
from opmodel import Machine

def main():
    cfg=json.load(sys.stdin);started=time.monotonic()
    c48,h48=leaders(P48);c32,h32=leaders(P32)
    l48=syndrome_luts(P48);l32=syndrome_luts(P32)
    s48=Fraction(sum(a*9**(48-w) for w,a in enumerate(h48)),10**48)
    s32=Fraction(sum(a*9**(32-w) for w,a in enumerate(h32)),10**32)
    success=s48**41*s32
    scale=10**100
    S=success.numerator*scale//success.denominator
    B=1-Fraction(S,scale)
    assert 1-success<=B<1
    manifest=json.loads((HERE/'submission.json').read_text())
    score=manifest['claims'][0]['score']
    assert B==Fraction(int(score['numerator']),int(score['denominator']))
    assert END < 1<<32 and W_ACTUAL==1804
    assert W_ACTUAL**1000<=2**10817
    assert hashlib.sha256((HERE/'opmodel.py').read_bytes()).hexdigest()=='9b440361c3fede634b8021705f6afd62e964c1d996045c91a16228acee296903'
    rng=random.Random(cfg['seed'])
    masks=[0,(1<<2000)-1]+[rng.getrandbits(2000) for _ in range(cfg['random_words'])]
    messages=[0,(1<<1000)-1]+[rng.getrandbits(1000) for _ in range(2)]
    masks += [encode(u) for u in messages]
    counts=[]
    for y in masks:
        m=Machine([(y>>i)&1 for i in range(2000)])
        m.mem=TableMemory(m.mem,c48,c32,l48,l32)
        decode(m); counts.append(m.ops)
        expected=0
        for b in range(42):
            k,n,start,dest,corr,lut=((24,48,48*b,24*b,c48,l48) if b<41 else (16,32,1968,984,c32,l32))
            word=(y>>start)&((1<<n)-1);s=0
            for j in range(n//8):s^=lut[256*j+((word>>(8*j))&255)]
            expected |= ((word^corr[s])&((1<<k)-1))<<dest
        actual=sum(bit<<i for i,bit in enumerate(m.output_bits()))
        assert actual==expected and m.ops==W_ACTUAL
    for u in messages:
        y=encode(u);m=Machine([(y>>i)&1 for i in range(2000)])
        m.mem=TableMemory(m.mem,c48,c32,l48,l32);decode(m)
        assert sum(bit<<i for i,bit in enumerate(m.output_bits()))==u
    # Exhaustively verify every byte syndrome-LUT entry against its linear map.
    for parity,lut in [(P48,l48),(P32,l32)]:
        cols=parity+[1<<i for i in range(len(parity))]
        for chunk in range(len(lut)//256):
            for value in range(256):
                s=0
                for bit in range(8):
                    if value>>bit&1:s^=cols[8*chunk+bit]
                assert lut[256*chunk+value]==s
    print(json.dumps({'histogram48':h48,'histogram32':h32,'cosets48':sum(h48),'cosets32':sum(h32),
        'component_success48':{'numerator':str(s48.numerator),'denominator':str(s48.denominator)},
        'component_success32':{'numerator':str(s32.numerator),'denominator':str(s32.denominator)},
        'success_floor_scaled_10_pow_100':str(S),'B_numerator':str(B.numerator),'B_denominator':str(B.denominator),
        'B_decimal_approx':float(B),'whole_success_approx':float(success),'exact_error_le_B':True,
        'reference_counter_sha256':hashlib.sha256((HERE/'opmodel.py').read_bytes()).hexdigest(),
        'program_primitive_count':W_ACTUAL,'claimed_W':1804,'W_pow_q_le_2_pow_p':True,'L_numerator':10817,'L_denominator':1000,
        'distinct_sample_counts':sorted(set(counts)),'sample_words':len(masks),'all_sample_outputs_correct':True,
        'registers':3,'memory_words_below_END':END,'table48_sha256':hashlib.sha256(c48.tobytes()).hexdigest(),
        'table32_sha256':hashlib.sha256(c32.tobytes()).hexdigest(),'elapsed_seconds':time.monotonic()-started}))
if __name__=='__main__':main()
