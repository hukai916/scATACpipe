"use strict"

// console.log("inside: ", _sectionAll);

// Solution:
  // step0: initialization to match with default page
  // step1: track all radio paramenter changes
    //1-1 set all other parameters display to 'none'
    //1-2 monitor the display of each section
    //1-3 change relevant other parameter to 'inline-block' for each displayed section
  // step2: handle doublet-removal-algorithm parameter to auto hide amulet-rmsk-bed and amulet-autosomes if not amulet.


let referenceGenomeSourceCurr = 'null';
let archrGenomeSourceCurr = 'null';
let isGenomeIndexBwaAvailableCurr = 'null';
let isGenomeIndexChromapAvailableCurr = 'null';
let isGenomeIndexCellrangerAvailableCurr = 'null';

// step0: initialization
let _allOtherParameterArr = [...document.querySelectorAll('.param-input')];
_allOtherParameterArr.forEach((item, i) => { item.style.display = 'none'; });

let _defaultParameterArr = [...document.querySelectorAll('.if-core, .if-fastq, .if-ensembl, .if-default')];
_defaultParameterArr.forEach((item, i) => { item.style.display = "inline-block"; });

// step1:
let _allRadioArr = [...document.querySelectorAll('input[type=radio]')];
_allRadioArr.forEach((item, i) => {
  item.addEventListener('change', function(e) {
    // 1-0: clear red borders
    [...document.querySelectorAll(".parameters-pipeline .param-input-value")].forEach((item, i) => {
      item.childNodes[1].style.border = '';
    });
    // 1-1
    _allOtherParameterArr.forEach((item, i) => { item.style.display = 'none'; });
    // .if-core should always on:
    [...document.querySelectorAll(".if-core")].forEach((item, i) => { item.style.display = "inline-block"; });

    // 1-2
      // input type section always on
    inputTypeRadioArr.forEach(function(elem) {
      if (elem.checked === true) {
        inputTypeCurr = elem.value;
      }
    })
    if (inputTypeCurr == "fragment") {
      [...document.querySelectorAll(".if-fragment")].forEach((item, i) => {
        item.style.display = "inline-block";
      })
    } else if (inputTypeCurr == "fastq") {
      [...document.querySelectorAll(".if-fastq")].forEach((item, i) => {
        item.style.display = "inline-block";
      })
    }

    // reference geome source section
    if (document.querySelector('div[id=div-reference-genome-source]').style.display == "block") {
      let referenceGenomeSourceRadioArr = [...document.querySelectorAll('input[type=radio][name="reference-genome-source"]')];
      referenceGenomeSourceRadioArr.forEach(function(elem) {
        if (elem.checked === true) {
          referenceGenomeSourceCurr = elem.value;
        }
      })
      if (referenceGenomeSourceCurr == 'ensembl') {
        [...document.querySelectorAll('.if-ensembl')].forEach((item, i) => { item.style.display = 'inline-block'; })
      } else if (referenceGenomeSourceCurr == 'ucsc') {
        [...document.querySelectorAll('.if-ucsc')].forEach((item, i) => { item.style.display = 'inline-block'; })
      } else if (referenceGenomeSourceCurr == 'custom-fasta') {
        [...document.querySelectorAll('.if-custom-fasta')].forEach((item, i) => { item.style.display = 'inline-block'; })
      }

    }

    // archr genome source section
    if (document.querySelector('div[id=div-archr-genome-source]').style.display == 'block') {
      let archrGenomeSourceRadioArr = [...document.querySelectorAll('input[type=radio][name="archr-genome-source"]')];
      archrGenomeSourceRadioArr.forEach(function(elem) {
        if (elem.checked === true) {
          archrGenomeSourceCurr = elem.value;
        }
      })

      if (archrGenomeSourceCurr == 'archr-defaults') {
        [...document.querySelectorAll('.if-archr-defaults')].forEach((item, i) => { item.style.display = "inline-block"; });
      } else if (archrGenomeSourceCurr == 'build-from-fasta') {
        [...document.querySelectorAll('.if-build-from-fasta')].forEach((item, i) => { item.style.display = "inline-block"; });
      } else if (archrGenomeSourceCurr == 'r-objects') {
        [...document.querySelectorAll('.if-r-objects')].forEach((item, i) => { item.style.display = "inline-block"; });
      }

    }

    // preprocessing strategy section
    if (document.querySelector('div[id=div-preprocessing-strategy]').style.display == 'block') {
      let preprocessingStrategyRadioArr = [...document.querySelectorAll('input[type=radio][name="preprocessing-strategy"]')];
      preprocessingStrategyRadioArr.forEach(function(elem) {
        if (elem.checked === true) {
          preprocessingStrategyCurr = elem.value;
        }
      })

      if (preprocessingStrategyCurr == 'default') {
        [...document.querySelectorAll('.if-default')].forEach((item, i) => { item.style.display = 'inline-block'; });
      } else if (preprocessingStrategyCurr == 'chromap') {
        [...document.querySelectorAll('.if-chromap')].forEach((item, i) => { item.style.display = 'inline-block'; });
      } else if (preprocessingStrategyCurr == 'default') {
        [...document.querySelectorAll('.if-10xgenomics')].forEach((item, i) => { item.style.display = 'inline-block'; });
      }
    }

    // is genome index bwa section
    if (document.querySelector('div[id=div-index-bwa]').style.display == 'block') {
      let indexBwaArr = [...document.querySelectorAll('[name="index-bwa-available"]')];
      indexBwaArr.forEach((item, i) => {
        if (item.checked === true) {
          isGenomeIndexBwaAvailableCurr = item.value;
        }
      })

      if (isGenomeIndexBwaAvailableCurr == 'index-bwa-yes') {
        [...document.querySelectorAll('.if-bwa-index')].forEach((item, i) => {
          item.style.display = 'inline-block';
        });
      }
    }

    // is genome index chromap section
    if (document.querySelector('div[id=div-index-chromap]').style.display == 'block') {
      let indexBwaArr = [...document.querySelectorAll('[name="index-chromap-available"]')];
      indexBwaArr.forEach((item, i) => {
        if (item.checked === true) {
          isGenomeIndexChromapAvailableCurr = item.value;
        }
      })

      if (isGenomeIndexChromapAvailableCurr == 'index-chromap-yes') {
        [...document.querySelectorAll('.if-chromap-index')].forEach((item, i) => {
          item.style.display = 'inline-block';
        });
      }
    }

    // is genome index cellranger section
    if (document.querySelector('div[id=div-index-cellranger]').style.display == 'block') {
      let indexBwaArr = [...document.querySelectorAll('[name="index-cellranger-available"]')];
      indexBwaArr.forEach((item, i) => {
        if (item.checked === true) {
          isGenomeIndexCellrangerAvailableCurr = item.value;
        }
      })

      if (isGenomeIndexCellrangerAvailableCurr == 'index-cellranger-yes') {
        [...document.querySelectorAll('.if-cellranger-index')].forEach((item, i) => {
          item.style.display = 'inline-block';
        });
      }
    }

  })
});

// step2:
  // initialize doublet-removal-algorithm to reflect default settings
  // since default uses archr, no need to display the other two params
  [...document.querySelectorAll('#param-amulet-rmsk-bed, #param-amulet-autosomes')].forEach((item, i) => {
    item.style.display = 'none';
  });

  // update upon change
document.querySelector('div[id=param-doublet-removal-algorithm').addEventListener('change', function(e) {
  if (document.querySelector('#doublet-removal-algorithm').value == 'amulet') {
    [...document.querySelectorAll('#param-amulet-rmsk-bed, #param-amulet-autosomes')].forEach((item, i) => {
      item.style.display = 'inline-block';
    });
  } else {
    [...document.querySelectorAll('#param-amulet-rmsk-bed, #param-amulet-autosomes')].forEach((item, i) => {
      item.style.display = 'none';
    });
  }
})
