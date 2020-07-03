package data

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strings"
)

// Refgene UCSC RefGene
type Refgene struct {
	Chrom      string   `json:"chrom"`
	Strand     byte     `json:"strand"`
	Gene       string   `json:"gene"`
	EntrezID   int      `json:"entrez_id"`
	Transcript string   `json:"translate"`
	ExonStart  int      `json:"exon_start"`
	ExonEnd    int      `json:"exon_end"`
	CdsStart   int      `json:"cds_start"`
	CdsEnd     int      `json:"cds_end"`
	ExonStarts []int    `json:"exon_starts"`
	ExonEnds   []int    `json:"exon_ends"`
	Regions    Regions  `json:"regions"`
	Streams    Regions  `json:"streams"`
	Tag        string   `json:"tag"`
	Mrna       Sequence `json:"mrna"`
	Cdna       Sequence `json:"cdna"`
	Protein    Sequence `json:"protein"`
}

// GetSn 获取RefGene编号
func (refgene Refgene) GetSn() string {
	return fmt.Sprintf("%s|%s:%d:%d", refgene.Transcript, refgene.Chrom, refgene.ExonStart, refgene.ExonEnd)
}

// IsCmpl RefGene 是否是完整转录本
func (refgene Refgene) IsCmpl() bool {
	return refgene.Tag == "cmpl"
}

// GetNumericalPosition 获取数值位置
func (refgene Refgene) GetNumericalPosition() (int, int) {
	order, _ := GetChromByName(refgene.Chrom)
	start := order*ChromMaxLen + refgene.Streams[0].Start
	end := order*ChromMaxLen + refgene.Streams[1].End
	return start, end
}

// SetRegions 设置区域元件
func (refgene *Refgene) SetRegions() {
	var regions Regions
	exonNum := len(refgene.ExonStarts)
	for i := 0; i < exonNum; i++ {
		if i > 0 {
			region := Region{
				Start: refgene.ExonEnds[i-1] + 1,
				End:   refgene.ExonStarts[i] - 1,
				Typo:  "intron",
			}
			regions = append(regions, region)
		}
		var exonOrder int
		if refgene.Strand == '+' {
			exonOrder = i + 1
		} else {
			exonOrder = exonNum - 1
		}
		start, end := refgene.ExonStarts[i], refgene.ExonEnds[i]
		if refgene.CdsStart > end || refgene.CdsEnd < start ||
			refgene.CdsStart <= start && end <= refgene.CdsEnd {
			var typo string
			if refgene.CdsStart > end {
				if refgene.Strand == '+' {
					typo = "utr5"
				} else {
					typo = "utr3"
				}
			} else if refgene.CdsEnd < start {
				if refgene.Strand == '-' {
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
			if start < refgene.CdsStart && refgene.CdsStart < end {
				if refgene.Strand == '+' {
					utrTypo1 = "utr5"
				} else {
					utrTypo2 = "utr3"
				}
				cdsStart = refgene.CdsStart
			}
			if start < refgene.CdsEnd && refgene.CdsEnd < end {
				if refgene.Strand == '+' {
					utrTypo2 = "utr3"
				} else {
					utrTypo2 = "utr5"
				}
				cdsEnd = refgene.CdsEnd
			}
			if utrTypo1 != "" {
				regions = append(regions, Region{
					Start:     start,
					End:       refgene.CdsStart - 1,
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
					Start:     refgene.CdsEnd + 1,
					End:       end,
					Typo:      utrTypo2,
					ExonOrder: exonOrder,
				})
			}
		}
	}
	sort.Sort(regions)
	refgene.Regions = regions
}

// SetUpDownStream 设置上下游区域
func (refgene *Refgene) SetUpDownStream(upDownStreamLen int) {
	stream1 := Region{
		Start: refgene.ExonStart - upDownStreamLen,
		End:   refgene.ExonStart - 1,
	}
	stream2 := Region{
		Start: refgene.ExonEnd + 1,
		End:   refgene.ExonEnd + upDownStreamLen,
	}
	if refgene.Strand == '+' {
		stream1.Typo = "upstream"
		stream2.Typo = "downstream"
	} else {
		stream1.Typo = "upstream"
		stream2.Typo = "downstream"
	}
	refgene.Streams = Regions{stream1, stream2}
}

// SetSequence 设置序列信息
func (refgene *Refgene) SetSequence(mrna Sequence) {
	if !mrna.IsEmpty() {
		refgene.Mrna = mrna
		if refgene.Tag != "unk" {
			for _, region := range refgene.Regions {
				if region.Typo == "cds" {
					seq := refgene.Mrna.GetSeq(region.Start-refgene.ExonStart, region.End-region.Start+1)
					refgene.Cdna.Join([]Sequence{seq})
				}
			}
			if refgene.Strand == '-' {
				refgene.Cdna.Reverse()
			}
		}
		if !refgene.Cdna.IsEmpty() {
			refgene.Protein = refgene.Cdna.Translate(refgene.Chrom == "MT")
			if refgene.Protein.IsCmpl() {
				refgene.Tag = "cmpl"
			} else {
				refgene.Tag = "incmpl"
				refgene.Protein.Clear()
			}
		}
	}
}

// NewRefgene 读取一行RefGene Line, 创建Refgene
func NewRefgene(refgeneLine string) (refgene Refgene, err error) {
	field := strings.Split(refgeneLine, "\t")
	// Position
	var posInts, startInts, endInts []int
	if posInts, err = Strs2Ints(field[4:8]); err != nil {
		return
	}
	if startInts, err = Strs2Ints(strings.Split(strings.Trim(field[9], ","), ",")); err != nil {
		return
	}
	for i := range startInts {
		startInts[i]++
	}
	if endInts, err = Strs2Ints(strings.Split(strings.Trim(field[9], ","), ",")); err != nil {
		return
	}
	refgene = Refgene{
		Chrom:      strings.Replace(field[2], "chr", "", 1),
		Transcript: field[1],
		Strand:     field[3][0],
		Gene:       field[12],
		Tag:        field[13],
		// EntrezID:   GetGeneEntrezID(refgene.Gene),
		ExonStart:  posInts[0] + 1,
		ExonEnd:    posInts[1],
		CdsStart:   posInts[2] + 1,
		CdsEnd:     posInts[3],
		ExonStarts: startInts,
		ExonEnds:   endInts,
	}
	refgene.SetRegions()
	refgene.SetUpDownStream(Config.Param.UpDownStream)
	return
}

// Refgenes Refgene 切片
type Refgenes []Refgene

// ToSnMap 转为sn编号为Key的Map集合
func (refgenes Refgenes) ToSnMap() (refgeneMap map[string]Refgene) {
	for _, refgene := range refgenes {
		refgeneMap[refgene.GetSn()] = refgene
	}
	return
}

// ToChromMap 转为chrom染色体为Key的Map集合
func (refgenes Refgenes) ToChromMap() (refgeneMap map[string]Refgenes) {
	for _, refgene := range refgenes {
		if rgs, ok := refgeneMap[refgene.Chrom]; ok {
			refgeneMap[refgene.Chrom] = append(rgs, refgene)
		} else {
			refgeneMap[refgene.Chrom] = Refgenes{refgene}
		}
	}
	return
}

// SetEntrezidAndSequence 向Refgenes中添加Entrez ID和Sequence信息
func (refgenes *Refgenes) SetEntrezidAndSequence(ncbiGene NcbiGene, mrna Fasta) {
	log.Printf("start set entrez id and sequence to refgenes")
	for i, refgene := range *refgenes {
		refgene.EntrezID = ncbiGene.GetEntrezID(refgene.Gene)
		if sequence, ok := mrna[refgene.GetSn()]; ok {
			refgene.SetSequence(sequence)
		}
		(*refgenes)[i] = refgene
	}
}

func (refgenes Refgenes) Len() int {
	return len(refgenes)
}

func (refgenes Refgenes) Less(i, j int) bool {
	starti, endi := refgenes[i].GetNumericalPosition()
	startj, endj := refgenes[j].GetNumericalPosition()
	if starti == startj {
		return endi < endj
	}
	return starti < startj
}

func (refgenes Refgenes) Swap(i, j int) {
	refgenes[i], refgenes[j] = refgenes[j], refgenes[i]
}

// ReadRefgeneFiles 读取Refgene文件
func ReadRefgeneFiles(refgeneFiles []string, refgenesChan chan Refgenes) {
	log.Printf("start read %s\n", strings.Join(refgeneFiles, ","))
	refgenes := make(Refgenes, 0)
	for _, refgeneFile := range refgeneFiles {
		var fp *os.File
		fp, err := os.Open(refgeneFile)
		if err != nil {
			log.Fatal(err)
		}
		defer fp.Close()
		reader := bufio.NewReader(fp)
		for {
			var line string
			line, err = reader.ReadString('\n')
			if err != nil {
				if err == io.EOF {
					break
				} else {
					log.Fatal(err)
				}
			}
			line = strings.TrimSpace(line)
			if len(line) == 0 || line[0] == '#' {
				continue
			}
			var refgene Refgene
			refgene, err = NewRefgene(line)
			if err != nil {
				log.Fatal(err)
			}
			if refgene.Chrom == "M" || len(refgene.Chrom) > 2 {
				continue
			}
			refgenes = append(refgenes, refgene)
		}
	}
	sort.Sort(refgenes)
	refgenesChan <- refgenes
}

// WriteMrnaFile 根据reference和refgeneMap输出mRNA序列Fasta文件
func WriteMrnaFile(mrnaFile string, refgenes Refgenes, reference Fasta) {
	log.Printf("start write %s\n", mrnaFile)
	fp, err := os.Create(mrnaFile)
	if err == nil {
		log.Fatal(err)
	}
	defer fp.Close()
	for _, refgene := range refgenes {
		if chromSeq, ok := reference[refgene.Chrom]; ok {
			mrna := chromSeq.GetSeq(refgene.ExonStart-1, refgene.ExonEnd-refgene.ExonStart+1)
			if _, err := fp.WriteString(fmt.Sprintf(">%s\n%s\n", refgene.GetSn(), mrna)); err != nil {
				log.Fatal(err)
			}
		}
	}
}
