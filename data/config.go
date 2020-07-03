package data

import (
	"io/ioutil"

	"gopkg.in/yaml.v2"
)

// ChromMaxLen 最大染色体长度
var ChromMaxLen int = 1000

// Config 配置
var Config struct {
	DBFile struct {
		Reference string `yaml:"reference"`
		NcbiGene  string `yaml:"ncbi_gene"`
		Refgene   string `yaml:"refgene"`
		EnsMt     string `yaml:"ens_mt"`
		Cds       string `yaml:"cds"`
		Exon      string `yaml:"exon"`
		Mrna      string `yaml:"mrna"`
		Refidx    string `yaml:"refidx"`
	} `yaml:"db_file"`
	Param struct {
		UpDownStream int `yaml:"up_down_stream"`
		RefidxStep   int `yaml:"refidx_step"`
		SplicingLen  int `yaml:"splicing_len"`
	} `yaml:"param"`
	Chrom []struct {
		Name   string `yaml:"name"`
		Length int    `yaml:"length"`
	} `yaml:"chrom"`
}

// ReadConfigYAML 读取YAML配置文件
func ReadConfigYAML(yamlFile string) error {
	buffer, err := ioutil.ReadFile(yamlFile)
	if err != nil {
		return err
	}
	err = yaml.Unmarshal(buffer, &Config)
	if err != nil {
		return err
	}
	maxLen := 0
	for _, chrom := range Config.Chrom {
		if maxLen < chrom.Length {
			maxLen = chrom.Length
		}
	}
	for {
		if ChromMaxLen/maxLen < 1 {
			break
		}
		ChromMaxLen *= 1000
	}
	return nil
}
