package data

import (
	"bytes"
	"log"
	"strconv"
	"strings"
)

// NcbiGene NCBI GENE INFO
type NcbiGene struct {
	Symbol    map[string]int
	Synonyms  map[string]int
	SymbolFna map[string]int
}

// GetEntrezID 通过Gene symbol获取Extrez ID
func (ng NcbiGene) GetEntrezID(symbol string) int {
	if entrezID, ok := ng.Symbol[symbol]; ok {
		return entrezID
	}
	if entrezID, ok := ng.SymbolFna[symbol]; ok {
		return entrezID
	}
	if entrezID, ok := ng.Synonyms[symbol]; ok {
		return entrezID
	}
	return -1
}

// ReadNCBIGeneInfo 读取NCBI GENE INFO文件
func ReadNCBIGeneInfo(ncbiGeneInfoFile string, ncbiGeneChan chan NcbiGene) {
	log.Printf("start read %s\n", ncbiGeneInfoFile)
	var ncbiGene NcbiGene
	var lines [][]byte
	lines, err := ReadFile(ncbiGeneInfoFile)
	if err != nil {
		log.Fatal(err)
	}
	for _, line := range lines {
		line = bytes.TrimSpace(line)
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		field := strings.Split(string(line), "\t")
		var entrezID int
		entrezID, err = strconv.Atoi(field[1])
		if err != nil {
			log.Fatal(err)
		}
		symbol, synonyms, symbolFna := field[2], strings.Split(field[4], "|"), field[10]
		if symbol != "-" && symbol != "." {
			ncbiGene.Symbol[symbol] = entrezID
		}
		if symbolFna != "-" && symbolFna != "." {
			ncbiGene.SymbolFna[symbolFna] = entrezID
		}
		for _, synonyms := range synonyms {
			if synonyms != "-" && synonyms != "." {
				ncbiGene.Synonyms[synonyms] = entrezID
			}
		}
	}
	ncbiGeneChan <- ncbiGene
}
