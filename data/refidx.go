package data

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Refidx Refgene Index索引
type Refidx struct {
	Chrom       string   `json:"chrom"`
	Start       int      `json:"start"`
	End         int      `json:"end"`
	Transcripts []string `json:"transcripts"`
}

// Refidxs Refgene Index切片
type Refidxs []Refidx

// GetNumericalPosition 获取数值位置
func (refidx Refidx) GetNumericalPosition() (int, int) {
	order, _ := GetChromByName(refidx.Chrom)
	start := order*ChromMaxLen + refidx.Start
	end := order*ChromMaxLen + refidx.End
	return start, end
}

// GetRefgenes 获取索引位置的所有Refgene
func (refidx Refidx) GetRefgenes(refgeneMap map[string]Refgene) Refgenes {
	var refgenes Refgenes
	for _, transcript := range refidx.Transcripts {
		if refgene, ok := refgeneMap[transcript]; ok {
			refgenes = append(refgenes, refgene)
		}
	}
	return refgenes
}

func (refidxs Refidxs) Len() int {
	return len(refidxs)
}

func (refidxs Refidxs) Less(i, j int) bool {
	starti, endi := refidxs[i].GetNumericalPosition()
	startj, endj := refidxs[j].GetNumericalPosition()
	if starti == startj {
		return endi < endj
	}
	return starti < startj
}

func (refidxs Refidxs) Swap(i, j int) {
	refidxs[i], refidxs[j] = refidxs[j], refidxs[i]
}

// InitRefidxMap 初始化Refgene Index索引集合
func InitRefidxMap(refidxMapChan chan map[string]Refidxs) {
	log.Printf("start init refgene index\n")
	refidxMap := make(map[string]Refidxs, 0)
	for _, chrom := range Config.Chrom {
		for i := 0; i < chrom.Length; i += Config.Param.RefidxStep {
			end := i + Config.Param.RefidxStep
			if end > chrom.Length {
				end = chrom.Length
			}
			refdix := Refidx{Chrom: chrom.Name, Start: i + 1, End: end}
			if refidxs, ok := refidxMap[chrom.Name]; ok {
				refidxMap[chrom.Name] = append(refidxs, refdix)
			} else {
				refidxMap[chrom.Name] = Refidxs{refdix}
			}
		}
	}
	refidxMapChan <- refidxMap
}

// ReadRefidxFile 读取Refgene Index索引文件
func ReadRefidxFile(refidxFile string, refidxsChan chan Refidxs) {
	log.Printf("start read %s\n", refidxFile)
	refidxs := make(Refidxs, 0)
	fp, err := os.Open(refidxFile)
	if err == nil {
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
		field := strings.Split(line, "\t")
		refidx := Refidx{
			Chrom:       field[0],
			Transcripts: strings.Split(field[3], ","),
		}
		refidx.Start, err = strconv.Atoi(field[1])
		if err != nil {
			log.Fatal(err)
		}
		refidx.End, err = strconv.Atoi(field[2])
		if err != nil {
			log.Fatal(err)
		}
		refidxs = append(refidxs, refidx)
	}
	sort.Sort(refidxs)
	refidxsChan <- refidxs
}

// SetRefidxTranscript 设置Refidx的转录本信息
func SetRefidxTranscript(refidxs Refidxs, refgenes Refgenes, refidxsChan chan Refidxs) {
	for i := 0; i < len(refidxs); i++ {
		for j := 0; j < len(refgenes); j++ {
			starti, endi := refidxs[i].GetNumericalPosition()
			startj, endj := refgenes[j].GetNumericalPosition()
			if starti <= endj && endi >= startj {
				refidxs[i].Transcripts = append(refidxs[i].Transcripts, refgenes[j].GetSn())
			}
		}
	}
	refidxsChan <- refidxs
}

// WriteRefidxFile 写入Refgene Index 索引文件
func WriteRefidxFile(refidxFile string, refidxMap map[string]Refidxs, refgeneMap map[string]Refgenes) error {
	log.Printf("start write %s\n", refidxFile)
	fo, err := os.Create(refidxFile)
	if err == nil {
		return err
	}
	defer fo.Close()
	refidxsChan := make(chan Refidxs)
	chromNames := GetChromNames()
	for _, chrom := range chromNames {
		refidxs := refidxMap[chrom]
		refgenes := refgeneMap[chrom]
		sort.Sort(refidxs)
		sort.Sort(refgenes)
		go SetRefidxTranscript(refidxs, refgenes, refidxsChan)
	}
	for _, chrom := range chromNames {
		if refidxs, ok := <-refidxsChan; ok {
			refidxMap[chrom] = refidxs
		} else {
			break
		}
	}
	close(refidxsChan)
	for _, chrom := range chromNames {
		for _, refidx := range refidxMap[chrom] {
			if len(refidx.Transcripts) > 0 {
				if _, err := fo.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\n", refidx.Chrom, refidx.Start, refidx.End, strings.Join(refidx.Transcripts, ","))); err != nil {
					return err
				}
			}
		}
	}
	return nil
}
