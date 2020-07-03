package cnv

import (
	"fmt"
	"grandanno/data"
	"log"
	"os"
	"strings"
)

// Annotation CNV注释
type Annotation struct {
	Gene       string `json:"gene"`
	EntrezID   int    `json:"entrez_id"`
	Transcript string `json:"transcript"`
	Region     string `json:"region"`
	Function   string `json:"function"`
	Exons      []int  `json:"exons"`
}

// AddExon 新增Exon信息
func (anno *Annotation) AddExon(exon int) {
	anno.Exons = append(anno.Exons, exon)
}

// GetCds 获取CDS信息
func (anno *Annotation) GetCds() string {
	switch len(anno.Exons) {
	case 0:
		return "."
	case 1:
		return fmt.Sprintf("exon.%d", anno.Exons[0])
	default:
		min, max := anno.Exons[0], anno.Exons[0]
		for _, exon := range anno.Exons {
			if min > exon {
				min = exon
			}
			if max < exon {
				max = exon
			}
		}
		return fmt.Sprintf("exon.%d_%d", min, max)
	}
}

// Annotations CNV注释切片
type Annotations []Annotation

// IsSpecial 是否包含重要变异
func (annos Annotations) IsSpecial() bool {
	for _, anno := range annos {
		if strings.Contains(anno.Region, "exon") {
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
func (annos *Annotations) AnnoStream(cnv Cnv, refgenes data.Refgenes) {
	for _, refgene := range refgenes {
		for _, region := range refgene.Streams {
			if cnv.GetVariant().Start <= region.End && cnv.GetVariant().End >= region.Start {
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
func (annos *Annotations) AnnoGene(cnv Cnv, refgenes data.Refgenes) {
	var cmplAnnos, incmplAnnos, unkAnnos Annotations
	variant := cnv.GetVariant()
	for _, refgene := range refgenes {
		if variant.End >= refgene.ExonStart && variant.Start <= refgene.ExonEnd {
			anno := Annotation{
				Gene:       refgene.Gene,
				EntrezID:   refgene.EntrezID,
				Transcript: refgene.Transcript,
			}
			if cnv.GetType() == "DEL" {
				anno.Function = "Deletion"
			} else {
				anno.Function = "Duplication"
			}
			if refgene.Tag == "unk" {
				anno.Region = "unkCDS"
				unkAnnos = append(unkAnnos, anno)
			} else {
				for _, region := range refgene.Regions {
					if variant.Start <= region.End && variant.End >= region.Start {
						anno.Region = region.Typo
						if region.Typo == "cds" && refgene.Tag == "cmpl" {
							anno.AddExon(region.ExonOrder)
						}
					}
				}
				if refgene.Tag == "cmpl" {
					cmplAnnos = append(cmplAnnos, anno)
				} else {
					anno.Region = "incmplCDS"
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
func RunAnnotation(cnvs Cnvs, refgeneMap map[string]data.Refgene, refidxs data.Refidxs, outJSONFile string) {
	fp, err := os.Create(outJSONFile)
	if err != nil {
		log.Fatal(err)
	}
	defer fp.Close()
	for i := 0; i < len(cnvs); i++ {
		cnvPos1, cnvPos2 := cnvs[i].GetVariant().GetNumericalPosition()
		for j := 0; j < len(refidxs); j++ {
			refPos1, refPos2 := refidxs[j].GetNumericalPosition()
			if cnvPos2 < refPos1 {
				break
			} else if cnvPos1 > refPos2 {
				continue
			} else {
				annos := make(Annotations, 0)
				if cnvPos2 < refPos1 {
					annos.AnnoIntergeic()
				} else {
					refgenes := refidxs[j].GetRefgenes(refgeneMap)
					annos.AnnoGene(cnvs[i], refgenes)
					if len(annos) == 0 {
						annos.AnnoStream(cnvs[i], refgenes)
					}
					if len(annos) == 0 {
						annos.AnnoIntergeic()
					}
				}
				json, err := data.ConvertToJSON(map[string]interface{}{"snv": cnvs[i], "annotations": annos})
				if err != nil {
					log.Fatal(err)
				}
				if _, err := fp.WriteString(json + "\n"); err != nil {
					log.Fatal(err)
				}
			}
		}
	}
}
