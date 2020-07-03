package data

import (
	"bufio"
	"bytes"
	"io"
	"log"
	"os"
	"strings"
)

// Base 碱基
type Base = byte

// Sequence 序列
type Sequence string

// Fasta fasta数据
type Fasta map[string]Sequence

// GetOne2Three 氨基酸简称(1字符)转全称(3字符)
func GetOne2Three(base Base) string {
	return CodonTbl.AA1to3[base]
}

// String 转原生字符串
func (seq Sequence) String() string {
	return string(seq)
}

// Reverse 翻转序列
func (seq *Sequence) Reverse() {
	var buffer bytes.Buffer
	for i := seq.GetLen() - 1; i >= 0; i-- {
		buffer.WriteByte(seq.GetChar(i))
	}
	*seq = Sequence(buffer.String())
}

// GetLen 获取序列长度
func (seq Sequence) GetLen() int {
	return len(seq)
}

// GetChar 获取自定位置的碱基
func (seq Sequence) GetChar(index int) Base {
	return seq[index]
}

// GetSeq 截取子序列
func (seq Sequence) GetSeq(index int, len int) Sequence {
	if len < 0 || index+len >= seq.GetLen() {
		return seq[index:]
	}
	return seq[index : index+len]
}

// GetIndex 获取自定碱基的位置（下标索引）
func (seq Sequence) GetIndex(base Base) int {
	return strings.IndexByte(seq.String(), base)
}

// IsStartswith 是否已指定碱基序列开头
func (seq Sequence) IsStartswith(prefix Sequence) bool {
	return strings.HasPrefix(seq.String(), prefix.String())
}

// IsEndswith 是否已指定碱基序列结尾
func (seq Sequence) IsEndswith(suffix Sequence) bool {
	return strings.HasSuffix(seq.String(), suffix.String())
}

// IsEqual 是否为同一序列
func (seq Sequence) IsEqual(seq2 Sequence) bool {
	return seq == seq2
}

// IsEmpty 序列是否为空
func (seq Sequence) IsEmpty() bool {
	return seq == ""
}

// RemoveOne 删除第一个子序列
func (seq *Sequence) RemoveOne(substr Sequence) {
	*seq = Sequence(strings.Replace(seq.String(), substr.String(), "", 1))
}

// RemoveAll 删除所有子序列
func (seq *Sequence) RemoveAll(substr Sequence) {
	*seq = Sequence(strings.Replace(seq.String(), substr.String(), "", -1))
}

// Clear 清空字符串
func (seq *Sequence) Clear() {
	*seq = ""
}

// Join 拼接多条字符串
func (seq *Sequence) Join(sequences []Sequence) {
	var buffer bytes.Buffer
	for _, sequence := range sequences {
		buffer.WriteString(sequence.String())
	}
	*seq = Sequence(buffer.String())
}

// GetSnpSequence 根据snp信息获取变异后的序列
func (seq Sequence) GetSnpSequence(pos int, alt Base) Sequence {
	var buffer bytes.Buffer
	buffer.WriteString(seq.GetSeq(0, pos-1).String())
	buffer.WriteByte(alt)
	buffer.WriteString(seq.GetSeq(pos, -1).String())
	return Sequence(buffer.String())
}

// GetInsSequence 根据insertion信息获取变异后的序列
func (seq Sequence) GetInsSequence(pos int, alt Sequence) Sequence {
	var buffer bytes.Buffer
	buffer.WriteString(seq.GetSeq(0, pos).String())
	buffer.WriteString(alt.String())
	buffer.WriteString(seq.GetSeq(pos, -1).String())
	return Sequence(buffer.String())
}

// GetDelSequence 根据deletion信息获取变异后的序列
func (seq Sequence) GetDelSequence(lenL int, lenR int) Sequence {
	var buffer bytes.Buffer
	buffer.WriteString(seq.GetSeq(0, lenL).String())
	buffer.WriteString(seq.GetSeq(seq.GetLen()-lenR, lenR).String())
	return Sequence(buffer.String())
}

// Translate 翻译为氨基酸序列
func (seq Sequence) Translate(isMt bool) Sequence {
	var buffer bytes.Buffer
	for i := 0; i < seq.GetLen(); i += 3 {
		codon := CodonTbl.Codon
		if isMt {
			codon = CodonTbl.CodonMt
		}
		if aa, ok := codon[seq.GetSeq(i, 3).String()]; ok {
			buffer.WriteByte(aa)
		}
	}
	return Sequence(buffer.String())
}

// IsCmpl Protein：是否为完整的氨基酸序列即存在终止密码子
func (seq Sequence) IsCmpl() bool {
	return strings.IndexByte(seq.String(), '*') != -1
}

// GetOne2Tree Protein：将氨基酸序列的每个Base由简称转为全称
func (seq Sequence) GetOne2Tree() string {
	var buffer bytes.Buffer
	for i := 0; i < seq.GetLen(); i++ {
		buffer.WriteString(GetOne2Three(seq[i]))
	}
	return buffer.String()
}

// ReadFastaFile 读取Fasta文件
func ReadFastaFile(fastaFile string, fastaChan chan Fasta) {
	log.Printf("start read %s\n", fastaFile)
	fasta := make(Fasta, 0)
	fp, err := os.Open(fastaFile)
	if err == nil {
		log.Fatal(err)
	}
	defer fp.Close()
	reader := bufio.NewReader(fp)
	var name, seq bytes.Buffer
	for {
		var line []byte
		line, err = reader.ReadBytes('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatal(err)
			}
		}
		line = bytes.TrimSpace(line)
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' {
			if name.Len() != 0 {
				fasta[strings.Split(name.String(), " ")[0]] = Sequence(seq.String())
			}
			name.Reset()
			seq.Reset()
			name.Write(line[1:])
		} else {
			seq.Write(line)
		}
	}
	if name.Len() != 0 {
		fasta[strings.Split(name.String(), " ")[0]] = Sequence(seq.String())
	}
	fastaChan <- fasta
}
