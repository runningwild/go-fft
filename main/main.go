package main

import (
  "fmt"
  "fft"
)

func main() {
  v,c := fft.Plan_sub(2*3*5*7*11, []int{2,3,5,7,11})
  fmt.Printf("%d: %v\n", c, v)
}
