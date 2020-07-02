package data

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"errors"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

// GetChromNames 获取染色体名称列表
func GetChromNames() (chromList []string) {
	for _, chrom := range Config.Chrom {
		chromList = append(chromList, chrom.Name)
	}
	return
}

// GetChromByName 通过名称获取染色体信息
func GetChromByName(name string) (order int, length int) {
	for index, chrom := range Config.Chrom {
		if chrom.Name == name {
			return index + 1, chrom.Length
		}
	}
	panic(errors.New("Not Found: " + name))
}

// GetGeneEntrezID 通过Gene symbol获取Extrez ID
func GetGeneEntrezID(symbol string) int {
	if entrezID, ok := NcbiGene.Symbol[symbol]; ok {
		return entrezID
	}
	if entrezID, ok := NcbiGene.SymbolFna[symbol]; ok {
		return entrezID
	}
	if entrezID, ok := NcbiGene.Synonyms[symbol]; ok {
		return entrezID
	}
	return -1
}

// ReadFile 读取文件全部内容
func ReadFile(file string) (lines [][]byte, err error) {
	fp, err := os.Open(file)
	if err != nil {
		return
	}
	defer fp.Close()
	var content []byte
	if strings.HasSuffix(strings.ToLower(file), ".gz") {
		var reader *gzip.Reader
		if reader, err = gzip.NewReader(fp); err != nil {
			return
		}
		if content, err = ioutil.ReadAll(reader); err != nil {
			return
		}
	} else {
		reader := bufio.NewReader(fp)
		if content, err = ioutil.ReadAll(reader); err != nil {
			return
		}
	}
	return bytes.Split(content, []byte{'\n'}), err
}

// Strs2Ints 批量字符串转整形
func Strs2Ints(strs []string) (ints []int, err error) {
	for _, s := range strs {
		i, err := strconv.Atoi(s)
		if err != nil {
			return []int{}, err
		}
		ints = append(ints, i)
	}
	return ints, nil
}
