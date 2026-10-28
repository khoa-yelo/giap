/* GIPPA report tables: click a header to sort, type in the box to filter (M9-14 behaviour 4).
   The rules are those of gippa/report/tables.py: numeric sort when every non-NA cell is a
   number, else string; ties keep their order; NA last; filter = case-insensitive substring. */
(function () {
  "use strict";
  var NA = "NA";
  var NUM = /^[+-]?((\d+\.?\d*|\.\d+)([eE][+-]?\d+)?|INF)$/;
  function num(c) {  /* Number("INF") is NaN: ADR-0006 writes infinities as INF / -INF */
    var m = /^([+-]?)INF$/.exec(c);
    return m ? (m[1] === "-" ? -Infinity : Infinity) : Number(c);
  }
  function order(cells, descending) {
    var idx = cells.map(function (c, i) { return i; });
    var keep = idx.filter(function (i) { return cells[i] !== NA; });
    var na = idx.filter(function (i) { return cells[i] === NA; });
    var numeric = keep.every(function (i) { return NUM.test(cells[i]); });
    var sign = descending ? -1 : 1;
    keep.sort(function (a, b) {
      var x = numeric ? num(cells[a]) : cells[a];
      var y = numeric ? num(cells[b]) : cells[b];
      return x < y ? -sign : x > y ? sign : 0;
    });
    return keep.concat(na);
  }
  document.querySelectorAll("table[data-table]").forEach(function (table) {
    var body = table.tBodies[0];
    table.querySelectorAll("thead th").forEach(function (th, col) {
      th.addEventListener("click", function () {
        var rows = Array.prototype.slice.call(body.rows);
        var desc = th.dataset.sort === "asc";
        th.dataset.sort = desc ? "desc" : "asc";
        var cells = rows.map(function (r) { return r.cells[col] ? r.cells[col].textContent : ""; });
        order(cells, desc).forEach(function (i) { body.appendChild(rows[i]); });
      });
    });
    var box = document.querySelector('input[data-filter="' + table.dataset.table + '"]');
    if (!box) { return; }
    box.addEventListener("input", function () {
      var q = box.value.toLowerCase();
      Array.prototype.forEach.call(body.rows, function (r) {
        r.hidden = q !== "" && !Array.prototype.some.call(r.cells, function (c) {
          return c.textContent.toLowerCase().indexOf(q) >= 0;
        });
      });
    });
  });
})();
