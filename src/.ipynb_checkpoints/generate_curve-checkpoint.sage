######################################################################
#  generate_prime_order_curve.sage
#  Usage (from Sage prompt):   E, G, candidates, sea_runs, gen_time = generate_prime_order_curve(256)
######################################################################

from sage.all import (
    ZZ, sqrt, randint, random_prime, Integer,
    GF, EllipticCurve, parallel, walltime
)

def safety_checks(E, p, n):
    """
    Mathematical sanity checks. n = E.cardinality()
    """
    # 2. Hasse bound …
    t = p + 1 - n
    assert abs(t) <= 2*sqrt(p),    "Violates Hasse bound!"
    # 3. Curve is non-singular by construction
    # 4. Provable primality of n
    assert n.is_prime(),           "Curve order is not provable prime!"


def generate_prime_order_curve(k, max_trials=10_000, log_every=50):
    """
    Fast prime-order curve generator over GF(p), p a random k-bit prime.
    Samples (a,b) at random. Logs every `log_every` SEA calls.
    Returns (E, G, candidates, sea_runs, gen_time) where
      - E is the resulting EllipticCurve
      - G is a generator point on E
      - candidates is the number of (a,b) pairs tested
      - sea_runs is the number of SEA executions performed
      - gen_time is the total time (in seconds) taken to generate the curve
    """
    sea_runs = 0
    candidates = 0
    t_start = walltime()

    for trial in range(max_trials):
        # pick a fresh random k-bit prime
        p = random_prime(2**k, lbound=2**(k-1))
        F = GF(p)

        while True:
            candidates += 1
            # uniform random coefficients
            a = randint(0, p-1)
            b = randint(0, p-1)

            # skip singular curves
            if (4*a**3 + 27*b**2) % p == 0:
                continue

            sea_runs += 1
            E = EllipticCurve(F, [a, b])
            n = E.cardinality()
            # log progress
            if sea_runs % log_every == 0:
                print(f"[SEA runs: {sea_runs}, candidates tested: {candidates}]")

            # quick composite filter
            if not n.is_prime(proof=False):
                continue

            # full safety checks
            safety_checks(E, p, n)

            # successful curve found
            print(f"✔ Found curve after SEA runs: {sea_runs}, candidates tested: {candidates}")
            G = E.random_point()
            t_end = walltime()
            gen_time = t_end - t_start
            return E, G, candidates, sea_runs, gen_time

    raise RuntimeError(f"No prime-order curve found in {max_trials} field trials")

######################################################################
# Quick demo / benchmark
######################################################################
if __name__ == "__main__":
    bits = 128
    print(f"Searching for a {bits}-bit prime-order curve …")
    E, G, candidates, sea_runs, gen_time = generate_prime_order_curve(bits)
    print(f"Generation complete in {gen_time:.2f} s")
    print(f"Candidates tested: {candidates}")
    print(f"SEA executions: {sea_runs}")
    print(f"p  = {E.base_field().order()}")
    print(f"a  = {E.a4()}")
    print(f"b  = {E.a6()}")
    print(f"#E = {E.cardinality()}")
    print(f"G  = {G.xy()}")
