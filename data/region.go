package data

// Region 区域元件：Exon，Intron，UpStream， DownStream
type Region struct {
	Start     int    `json:"start"`
	End       int    `json:"end"`
	Typo      string `json:"typo"`
	ExonOrder int    `json:"exon_order"`
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
