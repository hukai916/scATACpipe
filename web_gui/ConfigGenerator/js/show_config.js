let btnShowConfig = document.querySelector(".btn-show-script button");
let selectedInputMaster = "";

// function to clear all previously highlighted input borders if any.
let clearAllRedBorder = function() {
  [...document.querySelectorAll(".parameters-pipeline .param-input-value")].forEach((item, i) => {
    item.childNodes[1].style.border = '';
  });
}

// function to create download with pure JS without backend
let downloadConfig = function(filename, text) {
    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', filename);
    element.style.display = 'none';
    document.body.appendChild(element);
    element.click();
    document.body.removeChild(element);
}

btnShowConfig.addEventListener("click", function(e) {
  // alert(eventSubmit);
  // Clear previous codes if any:
  document.querySelector(".R-scripts textarea").value = "";
  // First toggle the Script div to block:
  document.querySelector(".R-scripts").style.display = "block";
  // Also clear all previous red highlighted borders if any:
  clearAllRedBorder();
  // Use setTimeout to let user know that the scripts have been updated.
  setTimeout(()=>{
    // Collect variable info and display R codes:
    let inputArr = [...document.querySelectorAll(".parameters-pipeline .param-input-value")];
    let inputVisible = [];
    inputArr.forEach((item, i) => {
      if (isvisible(item)) {
        inputVisible.push(item);
      }
    });

    // Warn user if any input is empty:
    let emptyInputArr = [];
    let rtextarea = document.querySelector(".R-scripts textarea");
    inputVisible.forEach(function(elem) {
      if (elem.childNodes[1].value == '') {
        emptyInputArr.push(elem.childNodes[1]);
      }
    })

    if (emptyInputArr.length > 0) {
      emptyInputArr.forEach(elem => elem.style.border = "0.2rem red solid");
      rtextarea.value = "Please fill in the hightlighted blanks!";
      rtextarea.style =  "font-size: 1.2rem; color: red; padding: 0.5rem;";
    } else {
      // display parameters input by user:
      inputVisible = inputVisible.filter(item => item.childNodes[1].name !== 'profile');

      let displayConfig = '// Save the following configurations into custom.config\n';
      displayConfig = displayConfig + "\tparams {\n"
      // figure out preprocessing strategy
      let _preprocessing = 'false';
      if (document.querySelector('div[id=div-preprocessing-strategy]').style.display == 'block') {
        let preprocessingStrategyRadioArr = [...document.querySelectorAll('input[type=radio][name="preprocessing-strategy"]')];
        preprocessingStrategyRadioArr.forEach(function(elem) {
          if (elem.checked === true) {
            preprocessingStrategyCurr = elem.value;
          }
        })
        _preprocessing = preprocessingStrategyCurr;
      }
      if (_preprocessing != 'false') {
        displayConfig = displayConfig + '\t\tpreprocess = ' + '\'' + _preprocessing + '\'\n';
      }

      let _tem = '';
      inputVisible.forEach((item, i) => {
        name = item.childNodes[1].name.replaceAll("-", "_");
        value = item.childNodes[1].value;
        value = value.replace(/(^"|"$|^'|'$)/g, ''); // get rid of leading/trailing quotes
        if (value != 'false' && value != 'true') {
          _tem = '\t\t' + name + ' = ' + '\'' + value + '\'';
        } else {
          _tem = _tem = '\t\t' + name + ' = ' + value;
        }
        displayConfig = displayConfig + _tem + "\n";
      });
      displayConfig = displayConfig + '\t}\n\n';
      displayConfig = displayConfig + "// Use this command: ";

      // figure out profile settings
      let _profileStr = '';
      let _profileArr = [];
      let _profile = [...document.querySelectorAll('input[name="profile"]')];
      _profile.forEach((item, i) => {
        if (item.checked) {
          _profileArr.push(item.value);
        }
      });

      if (_profileArr.length > 0) {
        _profileStr = _profileArr.join(',');
        displayConfig = displayConfig + " nextflow run main.nf -c custom.config -profile " + _profileStr;
      } else {
        displayConfig = displayConfig + " nextflow run main.nf -c custom.config";
      }

      rtextarea.style = "font-size: 1.2rem; color: black; padding: 0.5rem;";
      rtextarea.value = displayConfig;

    }

    let txtValue = document.querySelector(".R-scripts textarea").value;
    let configFile = '';

    if ((txtValue.length > 100) && (clickSubmit)) {
      downloadConfig("custom.config", txtValue);
    }
    clickSubmit = 0;
  }, 250);

})

// clear Show R scripts div box whenever any input/selection changes
let inputArr = [...document.querySelectorAll("select"), ...document.querySelectorAll("input"), ...document.querySelectorAll("textarea")];
inputArr.forEach((item, i) => {
  item.addEventListener("change", function(e) {
    document.querySelector(".R-scripts textarea").value = "";
    // document.querySelector(".R-scripts").style.display = "none";
  })
});

// add click event to download button
let eventFire = function(el, etype) {
  if (el.fireEvent) {
    el.fireEvent('on' + etype);
  } else {
    var evObj = document.createEvent('Events');
    evObj.initEvent(etype, true, false);
    el.dispatchEvent(evObj);
  }
}

// Add click event to submit button
let clickSubmit = 0;
let btn_submit = document.querySelector(".btn-submit-job button");
btn_submit.addEventListener("click", function(e) {
  clickSubmit = 1;
  eventFire(document.querySelector(".btn-show-script button"), "click");
})
