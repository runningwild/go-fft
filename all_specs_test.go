package fft_test

import (
  "gospec"
  "testing"
)


func TestAllSpecs(t *testing.T) {
  r := gospec.NewRunner()
  r.AddSpec(NaiveSpec)
  gospec.MainGoTest(r, t)
}

