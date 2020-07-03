package data

import (
	"fmt"
)

// Variant 变异
type Variant struct {
	Chrom string   `json:"chrom"`
	Start int      `json:"start"`
	End   int      `json:"end"`
	Ref   Sequence `json:"ref"`
	Alt   Sequence `json:"alt"`
}

// GetSn 获取变异编号
func (variant Variant) GetSn() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", variant.Chrom, variant.Start, variant.End, variant.Ref, variant.Alt)
}

// GetNumericalPosition 获取数值位置
func (variant Variant) GetNumericalPosition() (int, int) {
	order, _ := GetChromByName(variant.Chrom)
	start := order*ChromMaxLen + variant.Start
	end := order*ChromMaxLen + variant.End
	return start, end
}

// ConvertSnv 标准化变异信息
func (variant *Variant) ConvertSnv() {
	if variant.Chrom == "M" {
		variant.Chrom = "MT"
	}
	if !variant.Ref.IsEmpty() || !variant.Alt.IsEmpty() && !variant.Ref.IsEqual(variant.Alt) {
		if variant.Ref.IsStartswith(variant.Alt) || variant.Ref.IsEndswith(variant.Alt) {
			if variant.Ref.IsStartswith(variant.Alt) {
				variant.Start += variant.Alt.GetLen()

			}
			variant.Ref.RemoveOne(variant.Alt)
			variant.Alt.Clear()
		} else if variant.Alt.IsStartswith(variant.Ref) || variant.Alt.IsEndswith(variant.Ref) {
			if variant.Alt.IsStartswith(variant.Ref) {
				variant.Start += len(variant.Ref) - 1
			} else {
				variant.Start += len(variant.Ref) - len(variant.Alt)
			}
			variant.Alt.RemoveOne(variant.Ref)
			variant.Ref.Clear()
		} else {
			var refRev, altRev Sequence
			var subLen int
			refRev = variant.Ref
			altRev = variant.Alt
			refRev.Reverse()
			altRev.Reverse()
			for i, subLen := 0, 0; i < variant.Ref.GetLen() && i < variant.Alt.GetLen(); i++ {
				if refRev.GetChar(i) != altRev.GetChar(i) {
					break
				}
				subLen++
			}
			variant.Ref = variant.Ref.GetSeq(0, variant.Ref.GetLen()-subLen)
			variant.Alt = variant.Alt.GetSeq(0, variant.Alt.GetLen()-subLen)
			for i, subLen := 0, 0; i < variant.Ref.GetLen() && i < variant.Alt.GetLen(); i++ {
				if variant.Ref.GetChar(i) != variant.Alt.GetChar(i) {
					break
				}
				subLen++
			}
			variant.Ref = variant.Ref.GetSeq(subLen, -1)
			variant.Alt = variant.Alt.GetSeq(subLen, -1)
			if subLen > 0 && variant.Ref.IsEmpty() {
				variant.Start += subLen - 1
			} else {
				variant.Start += subLen
			}
		}
	}
	if variant.Ref.IsEmpty() {
		variant.End = variant.Start
		variant.Ref = "-"
	} else {
		variant.End = variant.Start + variant.Ref.GetLen() - 1
	}
	if variant.Alt.IsEmpty() {
		variant.Alt = "-"
	}
}
