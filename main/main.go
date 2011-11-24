package main

import (
  "fft"
  "fftw"
  "fmt"
  "math"
  "cmath"
)

func show(a []complex128) {
  for i,v := range a {
    re := real(v)
    im := real(v)
    if re < 1e-9 && re > -1e-9 { re = 0 }
    if im < 1e-9 && im > -1e-9 { im = 0 }
    fmt.Printf("%d:\t%2.4f\t%2.4f\n", i, re, im)
  }
  fmt.Printf("\n")
}

func twiddle(f, N int) complex128 {
  d := -2*math.Pi*fft.J * complex(float64(f), 0) / complex(float64(N),0)
  m := cmath.Exp(d)
  return m
}

func main2() {
  N := 16
  in := make([]complex128, N)
  for i := range in {
    in[i] = complex(math.Cos(float64(i) / float64(len(in)) * math.Pi * 2), 0)
  }
  out_fft := make([]complex128, N)
  out_fftw := make([]complex128, N)
  fftw.PlanDft1d(in, out_fftw, fftw.Forward, fftw.Estimate).Execute()
  show(in)
  show(out_fftw)
  fft.DFT(in, out_fft, 0, 2)
  fft.DFT(in, out_fft, 1, 2)

  for i := range in {
    in[i] = complex(0,0)
  }
  no2 := N / 2
  for i := 0; i < no2; i++ {
    in[i + 0] += out_fft[2*i + 0]
    in[i + 0] += out_fft[2*i + 1] * twiddle(i, N)
    in[i+no2] += out_fft[2*i + 0]
    in[i+no2] += out_fft[2*i + 1] * twiddle(i + no2, N)
  }
  show(in)
  fft.DFT(in[0:8], out_fft[0:8], 0, 1)
  fft.DFT(in[8: ], out_fft[8: ], 0, 1)
//  show(out_fft)
}

func main() {
  main7()
}

func main3() {
  N := 33
  in := make([]complex128, N)
  for i := range in {
    in[i] = complex(float64(i) / 12, float64(10 - i) / 13)
//    in[i] = complex(math.Cos(float64(i) / float64(len(in)) * math.Pi * 2), 0)
  }
  out_fft := make([]complex128, N)
  out_fftw := make([]complex128, N)
  fftw.PlanDft1d(in, out_fftw, fftw.Forward, fftw.Estimate).Execute()
  show(in)
  show(out_fftw)
  fft.DFT(in, out_fft, 0, 3)
  fft.DFT(in, out_fft, 1, 3)
  fft.DFT(in, out_fft, 2, 3)

  for i := range in {
    in[i] = complex(0,0)
  }
  no3 := N / 3
  for i := 0; i < no3; i++ {
    in[i + 0*no3] += out_fft[3*i + 0]
    in[i + 0*no3] += out_fft[3*i + 1] * twiddle(i, N)
    in[i + 0*no3] += out_fft[3*i + 2] * twiddle(2*i, N)

    in[i + 1*no3] += out_fft[3*i + 0]
    in[i + 1*no3] += out_fft[3*i + 1] * twiddle(i + no3, N)
    in[i + 1*no3] += out_fft[3*i + 2] * twiddle(2*(i + no3), N)

    in[i + 2*no3] += out_fft[3*i + 0]
    in[i + 2*no3] += out_fft[3*i + 1] * twiddle(i + 2*no3, N)
    in[i + 2*no3] += out_fft[3*i + 2] * twiddle(2*(i + 2*no3), N)

    // in[i + 0] += out_fft[2*i + 0]
    // in[i + 0] += out_fft[2*i + 1] * twiddle(i, N)
    // in[i+no2] += out_fft[2*i + 0]
    // in[i+no2] += out_fft[2*i + 1] * twiddle(i + no2, N)
  }
  show(in)
  fft.DFT(in[0:3], out_fft[0:3], 0, 1)
  fft.DFT(in[3:6], out_fft[3:6], 0, 1)
  fft.DFT(in[6:9], out_fft[6:9], 0, 1)
  fft.DFT(in[9: ], out_fft[9: ], 0, 1)
//  show(out_fft)
//  fft.DFT(in[0:8], out_fft[0:8], 0, 1)
//  fft.DFT(in[8: ], out_fft[8: ], 0, 1)
//  show(out_fft)
}

func butterfly(in,out []complex128, p,q int) {
  // len(v) == p*q
  for i := 0; i < p; i++ {
    for j := 0; j < q; j++ {
      for k := 0; k < q; k++ {
        out[i + j*p] += in[i*q + k] * twiddle(k * (i + j*p), len(in))
      }
    }
  }
  // no7 := N / 7
  // for i := 0; i < no7; i++ {
  //       // j                       k            k      j
  //   in[i + 0*no7] += out_fft[7*i + 0] * twiddle(0*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 1] * twiddle(1*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 2] * twiddle(2*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 3] * twiddle(3*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 4] * twiddle(4*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 5] * twiddle(5*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 6] * twiddle(6*(i + 0*no7), N)

  //   in[i + 1*no7] += out_fft[7*i + 0] * twiddle(0*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 1] * twiddle(1*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 2] * twiddle(2*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 3] * twiddle(3*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 4] * twiddle(4*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 5] * twiddle(5*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 6] * twiddle(6*(i + 1*no7), N)
}

func main7() {
  N := 35
  in := make([]complex128, N)
  for i := range in {
//    in[i] = complex(float64(i), float64(10 - i))
    in[i] = complex(math.Cos(float64(i) / float64(7+len(in)) * math.Pi * 2), 0)
  }
  out_fft := make([]complex128, N)
  out_fftw := make([]complex128, N)
  fftw.PlanDft1d(in, out_fftw, fftw.Forward, fftw.Estimate).Execute()
  show(in)
  show(out_fftw)
  fft.DFT(in, out_fft, 0, 7)
  fft.DFT(in, out_fft, 1, 7)
  fft.DFT(in, out_fft, 2, 7)
  fft.DFT(in, out_fft, 3, 7)
  fft.DFT(in, out_fft, 4, 7)
  fft.DFT(in, out_fft, 5, 7)
  fft.DFT(in, out_fft, 6, 7)

  for i := range in {
    in[i] = complex(0,0)
  }
  butterfly(out_fft, in, N/7, 7)
  // no7 := N / 7
  // for i := 0; i < no7; i++ {
  //   in[i + 0*no7] += out_fft[7*i + 0] * twiddle(0*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 1] * twiddle(1*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 2] * twiddle(2*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 3] * twiddle(3*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 4] * twiddle(4*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 5] * twiddle(5*(i + 0*no7), N)
  //   in[i + 0*no7] += out_fft[7*i + 6] * twiddle(6*(i + 0*no7), N)

  //   in[i + 1*no7] += out_fft[7*i + 0] * twiddle(0*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 1] * twiddle(1*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 2] * twiddle(2*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 3] * twiddle(3*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 4] * twiddle(4*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 5] * twiddle(5*(i + 1*no7), N)
  //   in[i + 1*no7] += out_fft[7*i + 6] * twiddle(6*(i + 1*no7), N)

  //   in[i + 2*no7] += out_fft[7*i + 0] * twiddle(0*(i + 2*no7), N)
  //   in[i + 2*no7] += out_fft[7*i + 1] * twiddle(1*(i + 2*no7), N)
  //   in[i + 2*no7] += out_fft[7*i + 2] * twiddle(2*(i + 2*no7), N)
  //   in[i + 2*no7] += out_fft[7*i + 3] * twiddle(3*(i + 2*no7), N)
  //   in[i + 2*no7] += out_fft[7*i + 4] * twiddle(4*(i + 2*no7), N)
  //   in[i + 2*no7] += out_fft[7*i + 5] * twiddle(5*(i + 2*no7), N)
  //   in[i + 2*no7] += out_fft[7*i + 6] * twiddle(6*(i + 2*no7), N)

  //   in[i + 3*no7] += out_fft[7*i + 0] * twiddle(0*(i + 3*no7), N)
  //   in[i + 3*no7] += out_fft[7*i + 1] * twiddle(1*(i + 3*no7), N)
  //   in[i + 3*no7] += out_fft[7*i + 2] * twiddle(2*(i + 3*no7), N)
  //   in[i + 3*no7] += out_fft[7*i + 3] * twiddle(3*(i + 3*no7), N)
  //   in[i + 3*no7] += out_fft[7*i + 4] * twiddle(4*(i + 3*no7), N)
  //   in[i + 3*no7] += out_fft[7*i + 5] * twiddle(5*(i + 3*no7), N)
  //   in[i + 3*no7] += out_fft[7*i + 6] * twiddle(6*(i + 3*no7), N)

  //   in[i + 4*no7] += out_fft[7*i + 0] * twiddle(0*(i + 4*no7), N)
  //   in[i + 4*no7] += out_fft[7*i + 1] * twiddle(1*(i + 4*no7), N)
  //   in[i + 4*no7] += out_fft[7*i + 2] * twiddle(2*(i + 4*no7), N)
  //   in[i + 4*no7] += out_fft[7*i + 3] * twiddle(3*(i + 4*no7), N)
  //   in[i + 4*no7] += out_fft[7*i + 4] * twiddle(4*(i + 4*no7), N)
  //   in[i + 4*no7] += out_fft[7*i + 5] * twiddle(5*(i + 4*no7), N)
  //   in[i + 4*no7] += out_fft[7*i + 6] * twiddle(6*(i + 4*no7), N)

  //   in[i + 5*no7] += out_fft[7*i + 0] * twiddle(0*(i + 5*no7), N)
  //   in[i + 5*no7] += out_fft[7*i + 1] * twiddle(1*(i + 5*no7), N)
  //   in[i + 5*no7] += out_fft[7*i + 2] * twiddle(2*(i + 5*no7), N)
  //   in[i + 5*no7] += out_fft[7*i + 3] * twiddle(3*(i + 5*no7), N)
  //   in[i + 5*no7] += out_fft[7*i + 4] * twiddle(4*(i + 5*no7), N)
  //   in[i + 5*no7] += out_fft[7*i + 5] * twiddle(5*(i + 5*no7), N)
  //   in[i + 5*no7] += out_fft[7*i + 6] * twiddle(6*(i + 5*no7), N)

  //   in[i + 6*no7] += out_fft[7*i + 0] * twiddle(0*(i + 6*no7), N)
  //   in[i + 6*no7] += out_fft[7*i + 1] * twiddle(1*(i + 6*no7), N)
  //   in[i + 6*no7] += out_fft[7*i + 2] * twiddle(2*(i + 6*no7), N)
  //   in[i + 6*no7] += out_fft[7*i + 3] * twiddle(3*(i + 6*no7), N)
  //   in[i + 6*no7] += out_fft[7*i + 4] * twiddle(4*(i + 6*no7), N)
  //   in[i + 6*no7] += out_fft[7*i + 5] * twiddle(5*(i + 6*no7), N)
  //   in[i + 6*no7] += out_fft[7*i + 6] * twiddle(6*(i + 6*no7), N)
  // }
  show(in)
  // fft.DFT(in[ 0: 7], out_fft[ 0: 7], 0, 1)
  // fft.DFT(in[ 7:14], out_fft[ 7:14], 0, 1)
  // fft.DFT(in[14:21], out_fft[14:21], 0, 1)
  // fft.DFT(in[21:28], out_fft[21:28], 0, 1)
  // fft.DFT(in[28:  ], out_fft[28:  ], 0, 1)
  // show(out_fft)
}
