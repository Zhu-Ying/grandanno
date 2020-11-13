package main

import (
	"grandanno/cnv"
	"grandanno/data"
	"grandanno/snv"
	"path"

	"github.com/spf13/cobra"
)

// Param 参数
var Param struct {
	Input          string
	Ouput          string
	Config         string
	DBPath         string
	SplicingLength int
}

// CorbaCMD 命令行参数解析
var CorbaCMD *cobra.Command

// prepareCMD 数据文件预处理
func prepareCMD() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "pre",
		Short: "预处理",
		Long:  "对数据库文件进行预处理得到注释所需文件",
		Run: func(cmd *cobra.Command, args []string) {
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
		},
	}
	cmd.Flags().StringVarP(&Param.DBPath, "config", "c", "config.yml", "配置文件")
	cmd.Flags().StringVarP(&Param.DBPath, "db_path", "d", "humandb", "数据库文件目录")
	return cmd
}

// annoGATKSNVCMD 注释GATK4 Call SNV的VCF结果文件
func annoGATKSNVCMD() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "SNV注释(GATK4)",
		Long:  "注释GATK4 Call SNV的VCF结果文件",
		Run: func(cmd *cobra.Command, args []string) {
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
		},
	}
	cmd.Flags().StringVarP(&Param.Config, "config", "c", "config.yml", "配置文件")
	cmd.Flags().StringVarP(&Param.DBPath, "db_path", "d", "humandb", "数据库文件目录")
	cmd.Flags().StringVarP(&Param.Input, "input", "i", "input.vcf", "输入vcf文件")
	cmd.Flags().StringVarP(&Param.Ouput, "output", "o", "output.json", "输出JSON文件")
	cmd.Flags().IntVarP(&Param.SplicingLength, "splicing_len", "s", -1, "预定义的剪接区域长度")
	return cmd
}

// annoXHMMCNVCMD 注释XHMM Call CNV的VCF结果文件
func annoXHMMCNVCMD() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cnv",
		Short: "CNV注释(XHMM)",
		Long:  "注释XHMM Call CNV的VCF结果文件",
		Run: func(cmd *cobra.Command, args []string) {
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
		},
	}
	cmd.Flags().StringVarP(&Param.Config, "config", "c", "config.yml", "配置文件")
	cmd.Flags().StringVarP(&Param.DBPath, "db_path", "d", "humandb", "数据库文件目录")
	cmd.Flags().StringVarP(&Param.Input, "input", "i", "input.vcf", "输入vcf文件")
	cmd.Flags().StringVarP(&Param.Ouput, "output", "o", "output.json", "输出JSON文件")
	return cmd
}

func init() {
	CorbaCMD = &cobra.Command{
		Use:   "grandanno",
		Short: "注释",
		Long:  "变异注释软件",
	}
	CorbaCMD.AddCommand(prepareCMD(), annoGATKSNVCMD(), annoXHMMCNVCMD())
	data.ReadConfigYAML(Param.Config)
	if Param.SplicingLength > 0 {
		data.Config.Param.SplicingLen = Param.SplicingLength
	}
}

func main() {
	CorbaCMD.Execute()
}
