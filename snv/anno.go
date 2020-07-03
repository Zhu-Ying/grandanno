package snv

import (
	"bytes"
	"grandanno/data"
	"log"
	"os"
	"strconv"
	"strings"
)

// Annotation 注释结果
type Annotation struct {
	Gene       string `json:"gene"`
	EntrezID   int    `json:"entrez_id"`
	Transcript string `json:"transcript"`
	Exon       string `json:"exon"`
	NaChange   string `json:"na_change"`
	AaChange   string `json:"aa_change"`
	Region     string `json:"region"`
	Function   string `json:"function"`
}

// SetExon 设置外显子信息
func (anno *Annotation) SetExon(exonOrder int) {
	var buffer bytes.Buffer
	buffer.WriteString("exon")
	buffer.WriteString(strconv.Itoa(exonOrder))
	anno.Exon = buffer.String()
}

// Annotations 注释结果切片
type Annotations []Annotation

// IsSpecial 判断注释结果是否为重要变异
func (annos Annotations) IsSpecial() bool {
	for _, anno := range annos {
		if strings.Contains(anno.Region, "splic") || strings.Contains(anno.Region, "exon") {
			return true
		}
	}
	return false
}

// AnnoIntergeic 注释基因间区
func (annos *Annotations) AnnoIntergeic() {
	*annos = append(*annos, Annotation{Region: "intergenic"})
}

// AnnoStream 注释上下游区
func (annos *Annotations) AnnoStream(snv Snv, refgenes data.Refgenes) {
	for _, refgene := range refgenes {
		for _, region := range refgene.Streams {
			variant := snv.GetVariant()
			if variant.Start <= region.End && variant.End >= region.Start {
				*annos = append(*annos, Annotation{
					Gene:       refgene.Gene,
					EntrezID:   refgene.EntrezID,
					Transcript: refgene.Transcript,
					Region:     region.Typo,
				})
				break
			}
		}

	}
}

// AnnoGene 注释基因区
func (annos *Annotations) AnnoGene(snv Snv, refgenes data.Refgenes, splicingLen int) {
	var cmplAnnos, incmplAnnos, unkAnnos Annotations
	for _, refgene := range refgenes {
		variant := snv.GetVariant()
		if variant.End >= refgene.ExonStart && variant.Start <= refgene.ExonEnd {
			anno := Annotation{
				Gene:       refgene.Gene,
				EntrezID:   refgene.EntrezID,
				Transcript: refgene.Transcript,
			}
			if refgene.Tag == "unk" {
				anno.Region = "unkCDS"
				unkAnnos = append(unkAnnos, anno)
			} else {
				switch snv.GetType() {
				case "del":
					anno.AnnoDel(snv, refgene, splicingLen)
				case "ins":
					anno.AnnoIns(snv, refgene, splicingLen)
				case "snp":
					anno.AnnoSnp(snv, refgene, splicingLen)
				default:
					anno.AnnoSnp(snv, refgene, splicingLen)
				}
				if refgene.IsCmpl() {
					cmplAnnos = append(cmplAnnos, anno)
				} else {
					anno.Function = "incmplCDS"
					incmplAnnos = append(incmplAnnos, anno)
				}
			}
		}
	}
	*annos = append(*annos, cmplAnnos...)
	if !(*annos).IsSpecial() {
		*annos = append(*annos, incmplAnnos...)
	}
	if len(*annos) == 0 {
		*annos = append(*annos, unkAnnos...)
	}
}

// RunAnnotation 运行注释
func RunAnnotation(snvs Snvs, refgeneMap map[string]data.Refgene, refidxs data.Refidxs, outJSONFile string) {
	log.Printf("start run annotation of snv\n")
	fp, err := os.Create(outJSONFile)
	if err != nil {
		log.Fatal(err)
	}
	defer fp.Close()
	for i, j := 0, 0; i < len(snvs) && j < len(refidxs); {
		snvPos1, snvPos2 := snvs[i].GetVariant().GetNumericalPosition()
		refPos1, refPos2 := refidxs[j].GetNumericalPosition()
		if snvPos1 > refPos2 {
			j++
		} else {
			annos := make(Annotations, 0)
			if snvPos2 < refPos1 {
				annos.AnnoIntergeic()
			} else {
				refgenes := refidxs[j].GetRefgenes(refgeneMap)
				annos.AnnoGene(snvs[i], refgenes, data.Config.Param.SplicingLen)
				if len(annos) == 0 {
					annos.AnnoStream(snvs[i], refgenes)
				}
				if len(annos) == 0 {
					annos.AnnoIntergeic()
				}
			}
			json, err := data.ConvertToJSON(map[string]interface{}{"snv": snvs[i], "annotations": annos})
			if err != nil {
				log.Fatal(err)
			}
			if _, err := fp.WriteString(json + "\n"); err != nil {
				log.Fatal(err)
			}
			i++
		}
	}
}
