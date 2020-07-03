package cnv

import "grandanno/data"

// Cnv CNV接口
type Cnv interface {
	GetVariant() data.Variant
	GetType() string
}

// Cnvs CNV接口切片
type Cnvs []Cnv

func (cnvs Cnvs) Len() int {
	return len(cnvs)
}

func (cnvs Cnvs) Less(i, j int) bool {
	starti, endi := cnvs[i].GetVariant().GetNumericalPosition()
	startj, endj := cnvs[j].GetVariant().GetNumericalPosition()
	if starti == startj {
		return endi < endj
	}
	return starti < startj
}

func (cnvs Cnvs) Swap(i, j int) {
	cnvs[i], cnvs[j] = cnvs[j], cnvs[i]
}
