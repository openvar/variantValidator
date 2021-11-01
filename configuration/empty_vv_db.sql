-- MySQL dump 10.13  Distrib 5.5.62, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: validator
-- ------------------------------------------------------
-- Server version	5.5.62-0ubuntu0.14.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `LRG_RSG_lookup`
--

DROP TABLE IF EXISTS `LRG_RSG_lookup`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `LRG_RSG_lookup` (
  `id` int(25) unsigned NOT NULL AUTO_INCREMENT,
  `lrgID` varchar(25) NOT NULL DEFAULT '',
  `hgncSymbol` varchar(25) NOT NULL,
  `RefSeqGeneID` varchar(25) NOT NULL DEFAULT '',
  `status` text NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `LRG_proteins`
--

DROP TABLE IF EXISTS `LRG_proteins`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `LRG_proteins` (
  `id` int(25) unsigned NOT NULL AUTO_INCREMENT,
  `LRGproteinID` varchar(25) DEFAULT NULL,
  `RefSeqProteinID` varchar(25) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `LRG_transcripts`
--

DROP TABLE IF EXISTS `LRG_transcripts`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `LRG_transcripts` (
  `id` int(25) unsigned NOT NULL AUTO_INCREMENT,
  `LRGtranscriptID` varchar(25) NOT NULL DEFAULT '',
  `RefSeqTranscriptID` varchar(25) NOT NULL DEFAULT '',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `refSeqGene_loci`
--

DROP TABLE IF EXISTS `refSeqGene_loci`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `refSeqGene_loci` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `refSeqGeneID` varchar(50) NOT NULL,
  `refSeqChromosomeID` varchar(500) NOT NULL,
  `genomeBuild` varchar(10) NOT NULL,
  `startPos` int(50) NOT NULL,
  `endPos` int(50) NOT NULL,
  `orientation` varchar(5) NOT NULL,
  `totalLength` int(50) NOT NULL,
  `chrPos` varchar(20) NOT NULL,
  `rsgPos` varchar(20) NOT NULL,
  `entrezID` int(20) NOT NULL,
  `hgncSymbol` varchar(20) NOT NULL,
  `updated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `transcript_info`
--

DROP TABLE IF EXISTS `transcript_info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `transcript_info` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `refSeqID` varchar(50) NOT NULL,
  `description` varchar(500) NOT NULL,
  `transcriptVariant` TEXT NOT NULL,
  `currentVersion` varchar(50) NOT NULL,
  `hgncSymbol` varchar(20) NOT NULL,
  `utaSymbol` varchar(20) NOT NULL,
  `updated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  KEY `refi` (`refSeqID`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `stableGeneIds`
--

DROP TABLE IF EXISTS `stableGeneIds`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `stableGeneIds` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `hgnc_symbol` varchar(25) DEFAULT NULL,
  `hgnc_id` varchar(25) DEFAULT NULL,
  `entrez_id` varchar(25) DEFAULT NULL,
  `ensembl_gene_id` varchar(100) DEFAULT NULL,
  `ucsc_id` varchar(100) DEFAULT NULL,
  `omim_id` varchar(100) DEFAULT NULL,
  `vega_id` varchar(100) DEFAULT NULL,
  `ccds_ids` varchar(1000) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;


--
-- Table structure for table `version`
--

DROP TABLE IF EXISTS `version`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `version` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `current_version` varchar(25) DEFAULT NULL,
  `updated` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;
