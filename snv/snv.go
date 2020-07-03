package snv

import "grandanno/data"

// Snv SNV接口
type Snv interface {
	GetVariant() data.Variant
	GetType() string
}

// Snvs SNV切片
type Snvs []Snv

func (snvs Snvs) Len() int {
	return len(snvs)
}

func (snvs Snvs) Less(i, j int) bool {
	starti, endi := snvs[i].GetVariant().GetNumericalPosition()
	startj, endj := snvs[j].GetVariant().GetNumericalPosition()
	if starti == startj {
		return endi < endj
	}
	return starti < startj
}

func (snvs Snvs) Swap(i, j int) {
	snvs[i], snvs[j] = snvs[j], snvs[i]
}
