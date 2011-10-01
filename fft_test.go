package fft_test

import (
  . "gospec"
  "gospec"
  "math"
  "cmath"
  "fft"
  "fftw"
  "fmt"
)

func FFTWSpec(c gospec.Context) {
  c.Specify("Check agains fftw", func() {
    N := 8 * 3 * 7 * 13
    in := make([]complex128, N)
    out_fft := make([]complex128, N)
    out_fftw := make([]complex128, N)
    for i := range in {
      in[i] = complex(float64(i / 1000 - 100), float64(i) / 10)
    }

    fftw.PlanDft1d(in, out_fftw, fftw.Forward, fftw.Estimate).Execute()
    fft.FFT(in, out_fft, 0, 1)

    for _,v := range out_fft {
      c.Expect(real(v), IsWithin(1e-9), real(v))
      c.Expect(imag(v), IsWithin(1e-9), imag(v))
    }
  })
}

func NaiveSpec(c gospec.Context) {
  in := make([]complex128, 8)
  for i := range in {
    in[i] = complex(math.Cos(float64(i) / float64(len(in)) * 2 * math.Pi), 0)
  }
  verify := func(out []complex128, name string) {
    c.Specify(name, func() {
//      fmt.Printf("name: %s\n", name)
//      for i := range out {
//        fmt.Printf("oot %d: %2.2f\n", i, out[i])
//      }
      for i := range out {
        mag,ang := cmath.Polar(out[i])
        if i == 1 || i == len(out) - 1 {
          c.Expect(mag, IsWithin(1e-9), float64(len(out))/2)
          if real(out[i]) < 0 {
            if ang < 0 {
              c.Expect(ang, IsWithin(1e-9), -math.Pi)
            } else {
              c.Expect(ang, IsWithin(1e-9), math.Pi)
            }
          } else {
            c.Expect(ang, IsWithin(1e-9), 0.0)
          }
        } else {
          c.Expect(mag, IsWithin(1e-9), 0.0)
        }
      }
    })
  }

  c.Specify("Test basic FFT", func() {
    in[3] = 3;
    out := make([]complex128, len(in))
    fft.FFT(in, out, 0, 1)
    for i := range out {
      fmt.Printf("oot %d: %2.2f\n", i, out[i])
    }
    fmt.Printf("\n\n")

    fft.DFT(in, out, 0, 1)
    for i := range out {
      fmt.Printf("oot %d: %2.2f\n", i, out[i])
    }
    fmt.Printf("\n\n")
  })

  c.Specify("Test basic DFT", func() {
    out := make([]complex128, len(in))

    fft.DFT(in, out, 0, 1)
    verify(out, "Start/Stride 0/1")

    in2 := make([]complex128, 2*len(in))
    for i := range in2 {
      in2[i] = in[i/2]
    }
    fft.DFT(in2, out, 0, 2)
    verify(out, "Start/Stride 0/2")


    fft.DFT(in2, out, 1, 2)
    verify(out, "Start/Stride 1/2")

    in5 := make([]complex128, len(in)*5)
    for i := range in5 {
      in5[i] = in[i/5]
    }
    for i := 0; i < 5; i++ {
      fft.DFT(in5, out, i, 5)
      verify(out, fmt.Sprintf("Start/Stride %d/%d", i, len(in5)/len(in)))
    }
  })
}