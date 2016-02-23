var DjangoResumable = function (options) {
    "use strict";
    var defaults, els;
    options = options || {};
    defaults = {
        csrfInputName: 'csrfmiddlewaretoken',
        urlAttribute: 'data-upload-url',
        progressDisplay: 'inline',
        errorListClass: 'errorlist',
        onFileError: this.onFileError,
        onFileAdded: this.onFileAdded,
        onFileSuccess: this.onFileSuccess,
        onFileProgress: this.onFileProgress,
        resumableOptions: {}
    };
    this.options = this.extend(defaults, options);
    this.csrfToken = document.querySelector('input[name=' + this.options.csrfInputName + ']').value;
    this.progressBars = {};
    els = document.querySelectorAll('input[' + this.options.urlAttribute + ']');
    this.each(els, function (el) {
        this.initField(el);
    });
};


DjangoResumable.prototype.each = function (elements, fn) {
    "use strict";
    var i, l;
    for (i = 0, l = elements.length; i < l; i += 1) {
        fn.apply(this, [elements[i]]);
    }
};

DjangoResumable.prototype.extend = function (target, source) {
    "use strict";
    var property;
    for (property in source) {
        if (source.hasOwnProperty(property)) {
            if (target[property] !== undefined && typeof source[property] === 'object') {
                target[property] = this.extend(target[property], source[property]);
            } else {
                target[property] = source[property];
            }

        }
    }
    return target;
};

DjangoResumable.prototype.getErrorList = function (el, create) {
    "use strict";
    var errorList = el.parentNode.previousSibling;
    while (errorList && errorList.tagName === undefined) {
        errorList = errorList.previousSibling;
    }
    if (errorList && !errorList.classList.contains(this.options.errorListClass)) {
        if (create === true) {
            errorList = document.createElement('ul');
            errorList.classList.add(this.options.errorListClass);
            el.parentNode.parentNode.insertBefore(errorList, el.parentNode);
        } else {
            errorList = null;
        }
    }
    return errorList;
};

DjangoResumable.prototype.initField = function (el) {
    "use strict";

    this.initResumable(el);

};


DjangoResumable.prototype.initProgressBar = function (file) {
    "use strict";

    // Create a bootstrap progress bar
    var progress = document.createElement('div');
    progress.className = 'progress';

    var bar = document.createElement('div');
    bar.className = "progress-bar primary_inv";
    bar.setAttribute('role', 'progressbar');
    bar.setAttribute('aria-valuenow', '0');
    bar.setAttribute('aria-valuemin', '0');
    bar.setAttribute('style', 'width:0%');
    bar.setAttribute('aria-valuemax', '100');
    progress.appendChild(bar);

    // Add it to the progressbar collection
    this.progressBars[file.fileName] = progress;

    // Add it to the DOM
    var el = file.container;

    el.parentNode.appendChild(progress);

    return progress;
};


DjangoResumable.prototype.initResumable = function (el) {
    "use strict";
    var elements = Array.prototype.slice.call(arguments),
        self = this,
        opts = {
            target: el.getAttribute(this.options.urlAttribute),
            query: {
                'csrfmiddlewaretoken': this.csrfToken
            }
        };

    opts = this.extend(this.options.resumableOptions, opts);
    var r = new Resumable(opts);
    r.assignBrowse(el);
    this.each(['fileAdded', 'fileProgress', 'fileSuccess', 'fileError'], function (eventType) {
        var callback = this.options['on' + eventType.substring(0, 1).toUpperCase() + eventType.substring(1)];
        r.on(eventType, function () {
            var args = arguments.length > 0 ? Array.prototype.slice.call(arguments) : [];
            callback.apply(self, [r].concat(args).concat(elements));
        });
    });

    return r;
};


DjangoResumable.prototype.onFileError = function (r, file, message, el) {
    "use strict";
    console.log(message);
    var errorList = this.getErrorList(el, true),
        error = document.createElement('li');
    error.innerHTML = message;
    if (errorList) {
        errorList.appendChild(error);
    }
};


DjangoResumable.prototype.onFileAdded = function (r, file, event, el) {
    "use strict";
    //var errorList = this.getErrorList(el);
    //if (errorList) {
    //    errorList.parentNode.removeChild(errorList);
    //}
    r.upload();

    // Add the progress bar
    this.initProgressBar(file);

    //Hide the input file
   el.style.display = 'none';


};


DjangoResumable.prototype.onFileSuccess = function (r, file, message, el) {
    "use strict";
    var progress = this.getProgressFromFile(file);

    $(progress).remove();
    delete this.progressBars[file.fileName];

    // If no more progress bars
    if ($.isEmptyObject(this.progressBars)) {

        $(el).show();
    }

};


DjangoResumable.prototype.onFileProgress = function (r, file, el) {
    "use strict";

    var progress = file.progress();
    var progressPctg = Math.floor(progress * 100);

    var progressBar = this.getProgressFromFile(file);

    var progressLine = progressBar.children[0];

    $(progressLine).css('width', progressPctg + '%').attr('aria-valuenow', progressPctg).html(file.fileName + ' - ' + progressPctg + '%');

};

DjangoResumable.prototype.getProgressFromFile = function (file) {

   return this.progressBars[file.fileName]

};