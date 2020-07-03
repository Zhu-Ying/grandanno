package cnv

import (
	"bytes"
	"grandanno/data"
	"log"
	"os"
	"strconv"
	"strings"
)

// XhmmCnv XHMM CNV 信息
type XhmmCnv struct {
	Variant     data.Variant `json:"Variant"`
	Information struct {
		MeanReadDepth         float64 `json:"mean_depth_depth"`
		MeanOriginalReadDepth float64 `json:"mean_original_depth"`
	} `json:"information"`
	OtherInfo []string `json:"other_info"`
}

// GetVariant 获取变异信息
func (cnv XhmmCnv) GetVariant() data.Variant {
	return cnv.Variant
}

// GetType 获取变异类型
func (cnv XhmmCnv) GetType() string {
	return strings.Trim(string(cnv.Variant.Alt), "<>")
}

// InitXhmmCnv 获取初始化XHMM CNV
func InitXhmmCnv(head []string, vcfLine string) (xhmmCnvMap map[string]Cnv, err error) {
	field := strings.Split(vcfLine, "\t")
	tmp1 := strings.Split(field[2], ":")
	chrom := tmp1[0] // important
	tmp2 := strings.Split(tmp1[1], "-")
	pos, err := data.Strs2Ints([]string{tmp2[0], tmp2[1]})
	if err != nil {
		return
	}
	start := pos[0]                      // important
	end := pos[1]                        // important
	ref := field[3]                      // important
	alts := strings.Split(field[4], ",") // important
	otherInfos := field[9:]              // important
	for i, otherInfo := range otherInfos {
		// var info Information
		sample := head[i]
		tmp3 := strings.Split(otherInfo, ":")
		lenght := len(tmp3)
		rd, err := data.Strs2Float([]string{tmp3[lenght-3], tmp3[lenght-1]})
		mrd := rd[0]  // important
		mord := rd[1] // important
		genotype, err := strconv.Atoi(tmp3[0])
		if err != nil {
			genotype = 0
		}
		if genotype == 0 {
			continue
		}
		xhmmCnv := XhmmCnv{
			Variant: data.Variant{
				Chrom: chrom,
				Start: start,
				End:   end,
				Ref:   data.Sequence(ref),
				Alt:   data.Sequence(alts[genotype-1]),
			},
			OtherInfo: otherInfos,
		}
		xhmmCnv.Information.MeanReadDepth = mrd
		xhmmCnv.Information.MeanOriginalReadDepth = mord
		xhmmCnvMap[sample] = xhmmCnv
	}
	return
}

// ReadXhmmVcfFile 读取XHMM VCF文件
func ReadXhmmVcfFile(vcfFile string, xhmmCnvMapChan chan map[string]Cnvs) {
	log.Printf("start read %s\n", vcfFile)
	xhmmCnvMap := make(map[string]Cnvs, 0)
	var head []string
	fp, err := os.Open(vcfFile)
	if err != nil {
		log.Fatal(err)
	}
	defer fp.Close()
	lines, err := data.ReadFile(vcfFile)
	if err != nil {
		log.Fatal(err)
	}
	for _, line := range lines {
		line = bytes.TrimSpace(line)
		if len(line) == 0 {
			continue
		}
		if line[0] == '#' {
			if bytes.HasPrefix(line, []byte("#CHROM")) {
				head = strings.Split(string(line), "\t")[9:]
			}
		} else {
			cnvMap, err := InitXhmmCnv(head, string(line))
			if err != nil {
				log.Fatal(err)
			}
			for sample, cnv := range cnvMap {
				if cnvs, ok := xhmmCnvMap[sample]; ok {
					xhmmCnvMap[sample] = append(cnvs, cnv)
				} else {
					xhmmCnvMap[sample] = Cnvs{cnv}
				}
			}
		}
	}
}
