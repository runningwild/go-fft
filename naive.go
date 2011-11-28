package fft

import (
  "math"
  "math/cmplx"
  "fmt"
)
func init() {
  fmt.Printf("")
}
const J = complex(0.0, 1.0)

// Performs a DFT or a partial DFT on in, storing the output in out
// if start == 0 && stride == 1 the DFT is done on the entire array,
// otherwise it is done on every stride elements, starting at start.
func DFT(in,out []complex128, start,stride int) {
  N := len(in) / stride
  if N == 1 {
    out[start] = in[start]
    return
  }
  factor := -2 * math.Pi * complex(0,1) / complex(float64(N),0)
  for k := 0; k < N; k++ {
    out[start + k*stride] = 0
    for n := start; n < len(in); n += stride {
      out[start + k*stride] += in[n] * cmplx.Exp(factor * complex(float64(k * (n / stride)), 0))
    }
  }
}

var n_factors []complex128
var kn_factors [][]complex128
const num_kn_factors = 64
const num_n_factors = num_kn_factors * num_kn_factors
func init() {
  n_factors = make([]complex128, num_n_factors)
  for i := range n_factors {
    n_factors[i] = -2 * math.Pi * J * complex(1.0 / float64(i), 0.0)
  }

  base_kn := make([]complex128, num_n_factors)
  kn_factors = make([][]complex128, num_kn_factors)
  for i := range kn_factors {
    kn_factors[i] = base_kn[i * num_kn_factors : (i+1) * num_kn_factors]
  }
  for i := range kn_factors {
    for j := range kn_factors {
      kfactor := n_factors[j] * complex(float64(i), 0.0)
      kn_factors[i][j] = cmplx.Exp(kfactor)
    }
  }
}

func twiddle(f, N int) complex128 {
  d := -2 * math.Pi * J * complex(float64(f), 0) / complex(float64(N),0)
  m := cmplx.Exp(d)
  return m
}

func butterfly(in,out []complex128, q, start, stride int) {
  // len(v) == p*q
  // TODO: this needs to be able to stride
  // start = 1, stride = 2
  //  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
  //     *     *     *     *     *     *     *     *     *     *
  // q = 5
  // p = 2
  // start + stride * (i + j*p)
  p := len(in) / q / stride
  for i := 0; i < p; i++ {
    for j := 0; j < q; j++ {
      out[start + stride * (i + j*p)] = 0
      for k := 0; k < q; k++ {
        out[start + stride * (i + j*p)] += in[start + stride * (i*q + k)] * twiddle((k * (i + j*p)), len(in) / stride)
      }
    }
  }
}

func factor(n int) []int {
  var f []int
  for i := 2; i*i <= n; i++ {
    for n%i == 0 {
      f = append(f, i)
      n /= i
      if n == 1 { return f }
    }
  }
  if n > 1 {
    f = append(f, n)
  }
  return f
}

type Plan struct {
  n int
  left,right *Plan
}

func (p *Plan) String() string {
  if p.left == nil && p.right == nil {
    return fmt.Sprintf("%d", p.n)
  }
  return fmt.Sprintf("%d {%v} {%v}", p.n, p.left, p.right)
}

// minimize pq(p+q)
// select some combination of factors and recurse
// returns min cost factors as a heap, min cost
func Plan_sub(n int, v []int) (*Plan,int64) {
  if len(v) == 1 {
    // TODO: decide whether we should really use v[0]*v[0] as the cost here, 
    // certainly not once we start using blustein
    return &Plan{ n : v[0]}, int64(v[0]*v[0])
  }
  if len(v) == 2 {
    return &Plan{ n, &Plan{ n : v[0] }, &Plan{ n : v[1] } }, int64(v[0]*v[1] * (v[0] + v[1]))
  }
  var a,b []int
  var best int64
  var best_plan *Plan
  for i := 1; i < 1<<uint(len(v)) - 1; i++ {
    a = a[0:0]
    b = b[0:0]
    var p int64 = 1
    for j := 0; j < len(v); j++ {
      if i & (1 << uint(j)) != 0 {
        a = append(a, v[j])
        p *= int64(v[j])
      } else {
        b = append(b, v[j])
      }
    }
    q := int64(n)/p
    f1s,cost1 := Plan_sub(int(p), a)
    f2s,cost2 := Plan_sub(int(q), b)
    total := q*p*(q+p) + cost1 + cost2

    if i == 1 || total < best {
      best = total
      best_plan = &Plan{ n, f1s, f2s }
    }
  }
  return best_plan, best
}


func fftHelper(in,out,temp []complex128, plan *Plan, start,stride int) {
  if plan.left == nil {
    DFT(in, out, start, stride)
    return
  }
  factor := plan.right.n
  for i := 0; i < factor; i++ {
//    sft_helper(in, out, factors, stride*factor)
    fftHelper(in, out, temp, plan.left, start + stride*i, stride*factor)
  }
  copy(temp, out)

  butterfly(temp, out, factor, start, stride)
}

func FFT(in,out []complex128) {
  temp := make([]complex128, len(in))
  plan,_ := Plan_sub(len(in), factor(len(in)))
  fftHelper(in, out, temp, plan, 0, 1)
}
