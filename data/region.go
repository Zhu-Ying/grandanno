package data

import "sort"

// Region 区域元件：Exon，Intron，UpStream， DownStream
type Region struct {
	Start     int
	End       int
	Typo      string
	ExonOrder int
}

// Regions 区域元件列表
type Regions []Region

// GetPrev 获取前一个区域元件
func (regions Regions) GetPrev(currentIndex int, strand byte) (Region, bool) {
	var index int
	if strand == '+' {
		index = currentIndex - 1
	} else {
		index = currentIndex + 1
	}
	if index >= 0 && index < len(regions) {
		return regions[index], true
	}
	return Region{}, false
}

// GetNext 获取后一个区域元件
func (regions Regions) GetNext(currentIndex int, strand byte) (Region, bool) {
	var index int
	if strand == '+' {
		index = currentIndex + 1
	} else {
		index = currentIndex - 1
	}
	if index >= 0 && index < len(regions) {
		return regions[index], true
	}
	return Region{}, false
}

func (regions Regions) Len() int {
	return len(regions)
}

func (regions Regions) Less(i, j int) bool {
	return regions[i].Start < regions[j].Start
}

func (regions Regions) Swap(i, j int) {
	regions[i], regions[j] = regions[j], regions[i]
}

// NewRegions 创建Regions
func NewRegions(position refgenePosition, strand byte) (regions Regions) {
	exonNum := len(position.ExonStarts)
	for i := 0; i < exonNum; i++ {
		if i > 0 {
			region := Region{
				Start: position.ExonEnds[i-1] + 1,
				End:   position.ExonStarts[i] - 1,
				Typo:  "intron",
			}
			regions = append(regions, region)
		}
		var exonOrder int
		if strand == '+' {
			exonOrder = i + 1
		} else {
			exonOrder = exonNum - 1
		}
		start, end := position.ExonStarts[i], position.ExonEnds[i]
		if position.CdsStart > end || position.CdsEnd < start ||
			position.CdsStart <= start && end <= position.CdsEnd {
			var typo string
			if position.CdsStart > end {
				if strand == '+' {
					typo = "utr5"
				} else {
					typo = "utr3"
				}
			} else if position.CdsEnd < start {
				if strand == '-' {
					typo = "utr5"
				} else {
					typo = "utr3"
				}
			} else {
				typo = "cds"
			}
			regions = append(regions, Region{
				Start:     start,
				End:       end,
				Typo:      typo,
				ExonOrder: exonOrder,
			})
		} else {
			utrTypo1, utrTypo2 := "", ""
			cdsStart, cdsEnd := start, end
			if start < position.CdsStart && position.CdsStart < end {
				if strand == '+' {
					utrTypo1 = "utr5"
				} else {
					utrTypo2 = "utr3"
				}
				cdsStart = position.CdsStart
			}
			if start < position.CdsEnd && position.CdsEnd < end {
				if strand == '+' {
					utrTypo2 = "utr3"
				} else {
					utrTypo2 = "utr5"
				}
				cdsEnd = position.CdsEnd
			}
			if utrTypo1 != "" {
				regions = append(regions, Region{
					Start:     start,
					End:       position.CdsStart - 1,
					Typo:      utrTypo1,
					ExonOrder: exonOrder,
				})
			}
			regions = append(regions, Region{
				Start:     cdsStart,
				End:       cdsEnd,
				Typo:      "cds",
				ExonOrder: exonOrder,
			})
			if utrTypo2 != "" {
				regions = append(regions, Region{
					Start:     position.CdsEnd + 1,
					End:       end,
					Typo:      utrTypo2,
					ExonOrder: exonOrder,
				})
			}
		}
	}
	sort.Sort(regions)
	return
}
