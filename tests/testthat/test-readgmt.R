library(BioQC)


context("Convert gmtlist (a memory-copy of a GMT file) to a list of signed gene sets")

testInputList <- list(list(name="GeneSetA_UP",genes=LETTERS[1:3]),
                      list(name="GeneSetA_DN", genes=LETTERS[4:6]),
                      list(name="GeneSetB", genes=LETTERS[2:4]),
                      list(name="GeneSetC_DN", genes=LETTERS[1:3]),
                      list(name="GeneSetD_UP", genes=LETTERS[1:3]))
outList.ignore <- gmtlist2signedGenesets(testInputList, nomatch="ignore")
outList.pos <- gmtlist2signedGenesets(testInputList, nomatch="pos")
outList.neg <- gmtlist2signedGenesets(testInputList, nomatch="neg")

exp.ignore <- list("GeneSetA"=list(pos=LETTERS[1:3], neg=LETTERS[4:6]),
                        "GeneSetB"=list(pos=NULL, neg=NULL),
                        "GeneSetC"=list(pos=NULL, neg=LETTERS[1:3]),
                        "GeneSetD"=list(pos=LETTERS[1:3], neg=NULL))

exp.pos <- list("GeneSetA"=list(pos=LETTERS[1:3], neg=LETTERS[4:6]),
                "GeneSetB"=list(pos=LETTERS[2:4], neg=NULL),
                "GeneSetC"=list(pos=NULL, neg=LETTERS[1:3]),
                "GeneSetD"=list(pos=LETTERS[1:3], neg=NULL))

exp.neg <- list("GeneSetA"=list(pos=LETTERS[1:3], neg=LETTERS[4:6]),
                "GeneSetB"=list(pos=NULL, neg=LETTERS[2:4]),
                "GeneSetC"=list(pos=NULL, neg=LETTERS[1:3]),
                "GeneSetD"=list(pos=LETTERS[1:3], neg=NULL))

test_that("gmtlist2signedGenesets, ignore non-matching genesets", {
              expect_equivalent(outList.ignore, exp.ignore)
          })

test_that("gmtlist2signedGenesets, ignore non-matching as positive", {
             expect_equivalent(outList.pos, exp.pos)
         })
test_that("gmtlist2signedGenesets, ignore non-matching as negative", {
             expect_equivalent(outList.neg, exp.neg)
         })
