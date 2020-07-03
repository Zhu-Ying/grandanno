package snv

import (
	"bytes"
	"grandanno/data"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

// GatkSnv GATK4 SNV 信息
type GatkSnv struct {
	Variant     data.Variant `json:"variant"`
	Information struct {
		Depth      int     `json:"depth"`
		Qual       float64 `json:"qual"`
		GatkFilter string  `json:"gatk_fitler"`
		Genotype   float64 `json:"genotype"`
		Ratio      float64 `json:"ratio"`
	} `json:"information"`
	OtherInfo string `json:"other_info"`
}

// GetVariant 获取变异信息
func (snv GatkSnv) GetVariant() data.Variant {
	return snv.Variant
}

// GetType 获取变异类型
func (snv GatkSnv) GetType() string {
	var typo string
	switch {
	case snv.Variant.Ref.GetChar(0) == '-':
		typo = "ins"
	case snv.Variant.Alt.GetChar(0) == '-':
		typo = "del"
	default:
		typo = "snp"
	}
	return typo
}

// InitGatkSnv 初始化GATK SNV
func InitGatkSnv(vcfLine string) (gatkSnvs Snvs, err error) {
	field := strings.Split(vcfLine, "\t")
	chrom := strings.Replace(field[0], "chr", "", 1)
	ref := field[3]
	alts := strings.Split(field[4], ",")
	gatkFilter := field[6]
	infoFeilds := strings.Split(field[7], ";")
	formatKeys := strings.Split(field[len(field)-2], ":")
	formatValues := strings.Split(field[len(field)-1], ":")
	qual, err := strconv.ParseFloat(field[5], 32)
	if err != nil {
		qual = -1
	}
	pos, err := strconv.Atoi(field[1])
	if err != nil {
		return
	}
	varCount := len(alts)
	// 获取基因型，覆盖深度、变异比率
	genotypes := make([]float64, varCount)
	ratios := make([]float64, varCount)
	depth := -1
	for i := 0; i < varCount; i++ {
		genotypes[i] = float64(-1)
		ratios[i] = float64(-1)
	}
	for _, info := range infoFeilds {
		if strings.HasPrefix(info, "AF=") {
			for i, gt := range strings.Split(info[3:], ",") {
				if i > varCount {
					break
				}
				if _gt, err := strconv.ParseFloat(gt, 32); err == nil {
					genotypes[i] = float64(_gt)
				}
			}
			break
		}
	}
	for i, key := range formatKeys {
		if key == "DP" {
			if dp, err := strconv.Atoi(formatValues[i]); err == nil {
				depth = dp
			}
		}
		if key == "AD" {
			sum := 0
			varCounts := make([]int, varCount)
			counts := strings.Split(formatValues[i], ",")
			for i := 0; i <= varCount && i < len(counts); i++ {
				if c, err := strconv.Atoi(counts[i]); err == nil {
					sum += c
					if i > 0 {
						varCounts[i-1] = c
					}
				}
			}
			if sum > 0 {
				for i := 0; i < varCount; i++ {
					ratios[i] = float64(varCounts[i]) / float64(sum)
				}
			}
		}
	}
	for i, alt := range alts {
		if alt == "*" {
			continue
		}
		gatkSnv := GatkSnv{
			Variant: data.Variant{
				Chrom: chrom,
				Start: pos,
				End:   0,
				Ref:   data.Sequence(ref),
				Alt:   data.Sequence(alt),
			},
			OtherInfo: vcfLine,
		}
		gatkSnv.Information.Depth = depth
		gatkSnv.Information.Qual = qual
		gatkSnv.Information.GatkFilter = gatkFilter
		gatkSnv.Information.Genotype = genotypes[i]
		gatkSnv.Information.Ratio = ratios[i]
		gatkSnv.Variant.ConvertSnv()
		gatkSnvs = append(gatkSnvs, gatkSnv)
	}
	return
}

// ReadGatkVcfFile 读取GATK VCF文件
func ReadGatkVcfFile(vcfFile string, gatkSnvsChan chan Snvs) {
	log.Printf("start read %s\n", vcfFile)
	gatkSnvs := make(Snvs, 0)
	fp, err := os.Open(vcfFile)
	if err != nil {
		log.Fatal(err)
	}
	defer fp.Close()
	// var lines [][]byte
	lines, err := data.ReadFile(vcfFile)
	if err != nil {
		log.Fatal(err)
	}
	for _, line := range lines {
		line = bytes.TrimSpace(line)
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		gatkSnv, err := InitGatkSnv(string(line))
		if err != nil {
			log.Fatal(err)
		}
		gatkSnvs = append(gatkSnvs, gatkSnv...)
	}
	sort.Sort(gatkSnvs)
	gatkSnvsChan <- gatkSnvs
}
