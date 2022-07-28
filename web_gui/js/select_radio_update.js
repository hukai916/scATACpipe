"use strict"

let inputTypeCurr = 'null';
// let referenceGenomeSourceCurr = 'null';
// let archrGenomeSourceCurr = 'null';
let preprocessingStrategyCurr = 'null';
// let isGenomeIndexBwaAvailableCurr = 'null';
// let isGenomeIndexChromapAvailableCurr = 'null';
// let isGenomeIndexCellrangerAvailableCurr = 'null';

// step0: initialization to match with the default settings in index.php
const _sectionAll = [...document.querySelectorAll('div[id=div-reference-genome-source], div[id=div-archr-genome-source], div[id=div-preprocessing-strategy], div[id=div-index-bwa], div[id=div-index-chromap], div[id=div-index-cellranger]')];
let _sectionToDisplay = [];
let inputTypeRadioArr = [...document.querySelectorAll('input[type=radio][name="input-type"]')];
inputTypeRadioArr.forEach(function(elem) {
  if (elem.checked === true) {
    inputTypeCurr = elem.value;
  }
})
if (inputTypeCurr == 'fragment') {
  _sectionToDisplay = [document.querySelector('div[id=div-archr-genome-source]')];
} else if (inputTypeCurr == 'fastq') {
  _sectionToDisplay = [...document.querySelectorAll('div[id=div-reference-genome-source], div[id=div-preprocessing-strategy], div[id=div-index-bwa]')];
}
_sectionAll.forEach((item, i) => { item.style.display = "none"; });
_sectionToDisplay.forEach((item, i) => { item.style.display = "block"; });

// step1: if input type changes
inputTypeRadioArr.forEach(function(elem){
  elem.addEventListener("change", function(e) {
    inputTypeRadioArr.forEach(function(elem) {
      if (elem.checked === true) {
        inputTypeCurr = elem.value;
      }
    })
    _sectionAll.forEach((item, i) => { item.style.display = "none"; });

    if (inputTypeCurr == 'fragment') {
      document.querySelector('div[id=div-archr-genome-source]').style.display = "block";
    } else if (inputTypeCurr == 'fastq') {
      // figure out which genome index section to display
      let preprocessingStrategyRadioArr = [...document.querySelectorAll('input[type=radio][name="preprocessing-strategy"]')];
      preprocessingStrategyRadioArr.forEach(function(elem) {
        if (elem.checked === true) {
          preprocessingStrategyCurr = elem.value;
        }
      })

      let _temIndex = '';
      if (preprocessingStrategyCurr == 'default') {
        _temIndex = 'div[id=div-index-bwa]';
      } else if (preprocessingStrategyCurr == 'chromap') {
        _temIndex = 'div[id=div-index-chromap]';
      } else if (preprocessingStrategyCurr == '10xgenomics') {
        _temIndex = 'div[id=div-index-cellranger]';
      }

      let _temString = 'div[id=div-reference-genome-source], div[id=div-preprocessing-strategy], ' + _temIndex;
      _sectionToDisplay = [...document.querySelectorAll(_temString)];

      _sectionAll.forEach((item, i) => { item.style.display = "none"; });
      _sectionToDisplay.forEach((item, i) => { item.style.display = "block"; });
    }
  })
})

// step2: if preprocessing strategy changes
  // need to update genone index section accordingly
let preprocessingStrategyRadioArr = [...document.querySelectorAll('input[type=radio][name="preprocessing-strategy"]')];
let _sectionGenomeIndexAll = [...document.querySelectorAll('div[id=div-index-bwa], div[id=div-index-chromap], div[id=div-index-cellranger]')]

preprocessingStrategyRadioArr.forEach(function(elem){
  elem.addEventListener('change', function(e) {
    preprocessingStrategyRadioArr.forEach(function(elem) {
      if (elem.checked === true) {
        preprocessingStrategyCurr = elem.value;
      }
    })
    _sectionGenomeIndexAll.forEach((item, i) => { item.style.display = 'none'; });

    console.log("preprocessing ", preprocessingStrategyCurr);

    if (preprocessingStrategyCurr == 'default') {
      document.querySelector('div[id=div-index-bwa]').style.display = 'block';
    } else if (preprocessingStrategyCurr == 'chromap') {
      document.querySelector('div[id=div-index-chromap]').style.display = 'block';
    } else if (preprocessingStrategyCurr == '10xgenomics') {
      document.querySelector('div[id=div-index-cellranger]').style.display = 'block';
    }
  })
})

// // step1: which sections to display
// let _sectionToDisplay = [];
//
//   // input type: always display
// let inputTypeRadioArr = [...document.querySelectorAll('input[type=radio][name="input-type"]')];
// inputTypeRadioArr.forEach(function(elem) {
//   if (elem.checked === true) {
//     inputTypeCurr = elem.value;
//   }
// })
//   // preprocessing strategy: used to determine which genome index section to display
// let preprocessingStrategyRadioArr = [...document.querySelectorAll('input[type=radio][name="preprocessing-strategy"]')];
// preprocessingStrategyRadioArr.forEach(function(elem) {
//   if (elem.checked === true) {
//     preprocessingStrategyCurr = elem.value;
//   }
// })
//
//
//
//   // figure out which sections to display
// if (inputTypeCurr == 'fragment') {
//   _sectionToDisplay = [document.querySelector('div[id=div-archr-genome-source]')];
// } else if (inputTypeCurr == 'fastq') {
//   _sectionToDisplay = [document.query]
// }
//
//
// // step2: which other parameter to display
