package data

import (
	"bytes"
	"strconv"
	"strings"
)

// NcbiGene NCBI GENE INFO
var NcbiGene struct {
	Symbol    map[string]int
	Synonyms  map[string]int
	SymbolFna map[string]int
}

// ReadNCBIGeneInfo 读取NCBI GENE INFO文件
func ReadNCBIGeneInfo(ncbiGeneInfoFile string) (err error) {
	NcbiGene.Symbol = make(map[string]int)
	NcbiGene.Synonyms = make(map[string]int)
	NcbiGene.SymbolFna = make(map[string]int)
	var lines [][]byte
	if lines, err = ReadFile(ncbiGeneInfoFile); err != nil {
		return
	}
	for _, line := range lines {
		line = bytes.TrimSpace(line)
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		field := strings.Split(string(line), "\t")
		var entrezID int
		if entrezID, err = strconv.Atoi(field[1]); err != nil {
			return
		}
		symbol, synonyms, symbolFna := field[2], strings.Split(field[4], "|"), field[10]
		if symbol != "-" && symbol != "." {
			NcbiGene.Symbol[symbol] = entrezID
		}
		if symbolFna != "-" && symbolFna != "." {
			NcbiGene.SymbolFna[symbolFna] = entrezID
		}
		for _, synonyms := range synonyms {
			if synonyms != "-" && synonyms != "." {
				NcbiGene.Synonyms[synonyms] = entrezID
			}
		}
	}
	return
}
