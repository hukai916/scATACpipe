// Below is what implemented in JQuery: $(object).is(':visible');
let isvisible = function (obj) {
  return obj.offsetWidth > 0 && obj.offsetHeight > 0;
}
