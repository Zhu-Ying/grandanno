package main

import (
	"flag"
	"grandanno/cnv"
	"grandanno/data"
	"grandanno/snv"
	"os"
	"path"
)

// Param 参数
var Param struct {
	Program        string
	Input          string
	Ouput          string
	Config         string
	DBPath         string
	SplicingLength int
	Help           bool
}

// runPrepare 数据文件预处理
func runPrepare() {
	referenceChan := make(chan data.Fasta, 0)
	go data.ReadFastaFile(path.Join(Param.DBPath, data.Config.DBFile.Reference), referenceChan)
	refgenesChan := make(chan data.Refgenes, 0)
	refgeneFiles := []string{path.Join(Param.DBPath, data.Config.DBFile.Refgene), path.Join(Param.DBPath, data.Config.DBFile.EnsMt)}
	go data.ReadRefgeneFiles(refgeneFiles, refgenesChan)
	refidxMapChan := make(chan map[string]data.Refidxs, 0)
	go data.InitRefidxMap(refidxMapChan)
	refgenes := <-refgenesChan
	data.WriteMrnaFile(path.Join(Param.DBPath, data.Config.DBFile.Mrna), refgenes, <-referenceChan)
	data.WriteRefidxFile(path.Join(Param.DBPath, data.Config.DBFile.Refidx), <-refidxMapChan, refgenes.ToChromMap())
}

// runPrepare 注释GATK SNV
func runAnnoGatkSnv() {
	ncbiGeneChan := make(chan data.NcbiGene, 0)
	go data.ReadNCBIGeneInfo(path.Join(Param.DBPath, data.Config.DBFile.NcbiGene), ncbiGeneChan)
	refgenesChan := make(chan data.Refgenes, 0)
	refgeneFiles := []string{path.Join(Param.DBPath, data.Config.DBFile.Refgene), path.Join(Param.DBPath, data.Config.DBFile.EnsMt)}
	go data.ReadRefgeneFiles(refgeneFiles, refgenesChan)
	mrnaChan := make(chan data.Fasta, 0)
	go data.ReadFastaFile(path.Join(Param.DBPath, data.Config.DBFile.Mrna), mrnaChan)
	refidxsChan := make(chan data.Refidxs)
	go data.ReadRefidxFile(path.Join(Param.DBPath, data.Config.DBFile.Refidx), refidxsChan)
	refgenes := <-refgenesChan
	refgenes.SetEntrezidAndSequence(<-ncbiGeneChan, <-mrnaChan)
	gatkSnvsChan := make(chan snv.Snvs)
	go snv.ReadGatkVcfFile(Param.Input, gatkSnvsChan)
	snv.RunAnnotation(<-gatkSnvsChan, refgenes.ToSnMap(), <-refidxsChan, Param.Ouput)
}

// runPrepare 注释XHMM CNV
func runAnnoXhmmCnv() {
	ncbiGeneChan := make(chan data.NcbiGene, 0)
	go data.ReadNCBIGeneInfo(path.Join(Param.DBPath, data.Config.DBFile.NcbiGene), ncbiGeneChan)
	refgenesChan := make(chan data.Refgenes, 0)
	refgeneFiles := []string{path.Join(Param.DBPath, data.Config.DBFile.Refgene), path.Join(Param.DBPath, data.Config.DBFile.EnsMt)}
	go data.ReadRefgeneFiles(refgeneFiles, refgenesChan)
	refidxsChan := make(chan data.Refidxs)
	go data.ReadRefidxFile(path.Join(Param.DBPath, data.Config.DBFile.Refidx), refidxsChan)
	xhmmCnvMapChan := make(chan map[string]cnv.Cnvs, 0)
	go cnv.ReadXhmmVcfFile(Param.Input, xhmmCnvMapChan)
	refgenes := <-refgenesChan
	refidxs := <-refidxsChan
	for sample, cnvs := range <-xhmmCnvMapChan {
		outJSONFile := path.Join(Param.Ouput + "." + sample + ".json")
		cnv.RunAnnotation(cnvs, refgenes.ToSnMap(), refidxs, outJSONFile)
	}
}

func init() {
	flag.StringVar(&Param.Program, "p", "gatk_snv", "The Program which is run:[prepare, gatk_snv, xhmm_cnv]")
	flag.StringVar(&Param.Input, "i", "", "Input File")
	flag.StringVar(&Param.Ouput, "o", "", "Output File")
	flag.StringVar(&Param.Config, "c", "", "Config YAML File")
	flag.StringVar(&Param.DBPath, "d", "./database", "Directory Path of Database")
	flag.IntVar(&Param.SplicingLength, "s", 15, "Length of Splicing Region")
	flag.BoolVar(&Param.Help, "h", false, "Help")
	flag.Parse()
	if Param.Help {
		flag.Usage()
		os.Exit(0)
	}
	if len(Param.Input) == 0 || len(Param.Ouput) == 0 || len(Param.Config) == 0 {
		flag.Usage()
		os.Exit(-1)
	}
	data.ReadConfigYAML(Param.Config)
	data.Config.Param.SplicingLen = Param.SplicingLength
}

func main() {
	switch Param.Program {
	case "prepare":
		runPrepare()
	case "gatk_snv":
		runAnnoGatkSnv()
	case "xhmm_snv":
		runAnnoXhmmCnv()
	default:
		flag.Usage()
	}
}
