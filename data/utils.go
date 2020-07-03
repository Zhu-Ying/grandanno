package data

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/json"
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
	ints = make([]int, len(strs))
	for i, s := range strs {
		if v, e := strconv.Atoi(s); e != nil {
			ints[i] = -1
			err = e
		} else {
			ints[i] = v
		}
	}
	return
}

// Strs2Float 批量字符串转浮点型
func Strs2Float(strs []string) (floats []float64, err error) {
	floats = make([]float64, len(strs))
	for i, s := range strs {
		if v, e := strconv.ParseFloat(s, 32); e != nil {
			floats[i] = float64(-1)
			err = e
		} else {
			floats[i] = v
		}
	}
	return
}

// ConvertToJSON 装换为JSON字符串
func ConvertToJSON(data interface{}) (string, error) {
	var buffer bytes.Buffer
	jsonEncoder := json.NewEncoder(&buffer)
	jsonEncoder.SetEscapeHTML(false)
	if err := jsonEncoder.Encode(data); err != nil {
		return "", err
	}
	return buffer.String(), nil
}
