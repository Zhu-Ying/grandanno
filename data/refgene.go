package data

import (
	"fmt"
	"strings"
)

type refgenePosition struct {
	ExonStart  int
	ExonEnd    int
	CdsStart   int
	CdsEnd     int
	ExonStarts []int
	ExonEnds   []int
}

// Refgene UCSC RefGene
type Refgene struct {
	Chrom      string
	Strand     byte
	Gene       string
	EntrezID   int
	Transcript string
	Position   refgenePosition
	Regions    Regions
	Streams    Regions
	Tag        string
	Mrna       Sequence
	Cdna       Sequence
	Protein    Sequence
}

// RefgeneMap Refgene Map集合
var RefgeneMap map[string]Refgene

// GetSn 获取RefGene编号
func (refgene Refgene) GetSn() string {
	return fmt.Sprintf("%s|%s:%d:%d", refgene.Transcript, refgene.Chrom, refgene.Position.ExonStart, refgene.Position.ExonEnd)
}

// IsCmpl RefGene 是否是完整转录本
func (refgene Refgene) IsCmpl() bool {
	return refgene.Tag == "cmpl"
}

// SetUpDownStream 设置上下游区域
func (refgene *Refgene) SetUpDownStream(upDownStreamLen int) {
	stream1 := Region{
		Start: refgene.Position.ExonStart - upDownStreamLen,
		End:   refgene.Position.ExonStart - 1,
	}
	stream2 := Region{
		Start: refgene.Position.ExonEnd + 1,
		End:   refgene.Position.ExonEnd + upDownStreamLen,
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
func (refgene *Refgene) SetSequence(sequence Sequence) {
	if !sequence.IsEmpty() {
		refgene.Mrna = sequence
		if refgene.Tag != "unk" {
			for _, region := range refgene.Regions {
				if region.Typo == "cds" {
					seq := refgene.Mrna.GetSeq(region.Start-refgene.Position.ExonStart, region.End-region.Start+1)
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
		EntrezID:   GetGeneEntrezID(refgene.Gene),
		Position: refgenePosition{
			ExonStart:  posInts[0] + 1,
			ExonEnd:    posInts[1],
			CdsStart:   posInts[2] + 1,
			CdsEnd:     posInts[3],
			ExonStarts: startInts,
			ExonEnds:   endInts,
		},
	}
	refgene.Regions = NewRegions(refgene.Position, refgene.Strand)
	return
}

func (refgeneDict RefgeneDict) SetUpDownStream(upDownStreamLen int) {
	for sn, refgene := range refgeneDict {
		refgene.SetUpDownStream(upDownStreamLen)
		refgeneDict[sn] = refgene
	}
}

func (refgeneDict RefgeneDict) SetSequence(mrna Fasta) {
	for sn, refgene := range refgeneDict {
		if sequence, ok := mrna[sn]; ok {
			refgene.SetSequence(sequence)
			refgeneDict[sn] = refgene
		}
	}
}
