(function($) {
	var num = function(value) {
		return parseInt(value, 10) || 0;
	};

	/**
	 * Sets or gets the values for min-width, min-height, max-width and
	 * max-height.
	 */
	$
			.each(
					[ 'min', 'max' ],
					function(i, name) {
						$.fn[name + 'Size'] = function(value) {
							var width, height;
							if (value) {
								if (value.width !== undefined) {
									this.css(name + '-width', value.width);
								}
								if (value.height !== undefined) {
									this.css(name + '-height', value.height);
								}
								return this;
							} else {
								width = this.css(name + '-width');
								height = this.css(name + '-height');
								// Apparently:
								// * Opera returns -1px instead of none
								// * IE6 returns undefined instead of none
								return {
									'width' : (name === 'max'
											&& (width === undefined
													|| width === 'none' || num(width) === -1) && Number.MAX_VALUE)
											|| num(width),
									'height' : (name === 'max'
											&& (height === undefined
													|| height === 'none' || num(height) === -1) && Number.MAX_VALUE)
											|| num(height)
								};
							}
						};
					});

	/**
	 * Returns whether or not an element is visible.
	 */
	$.fn.isVisible = function() {
		return this.is(':visible');
	};

	/**
	 * Sets or gets the values for border, margin and padding.
	 */
	$.each([ 'border', 'margin', 'padding' ],
			function(i, name) {
				$.fn[name] = function(value) {
					if (value) {
						if (value.top !== undefined) {
							this.css(name + '-top'
									+ (name === 'border' ? '-width' : ''),
									value.top);
						}
						if (value.bottom !== undefined) {
							this.css(name + '-bottom'
									+ (name === 'border' ? '-width' : ''),
									value.bottom);
						}
						if (value.left !== undefined) {
							this.css(name + '-left'
									+ (name === 'border' ? '-width' : ''),
									value.left);
						}
						if (value.right !== undefined) {
							this.css(name + '-right'
									+ (name === 'border' ? '-width' : ''),
									value.right);
						}
						return this;
					} else {
						return {
							top : num(this.css(name + '-top'
									+ (name === 'border' ? '-width' : ''))),
							bottom : num(this.css(name + '-bottom'
									+ (name === 'border' ? '-width' : ''))),
							left : num(this.css(name + '-left'
									+ (name === 'border' ? '-width' : ''))),
							right : num(this.css(name + '-right'
									+ (name === 'border' ? '-width' : '')))
						};
					}
				};
			});
})(jQuery);

(function() {
	jLayout = (typeof jLayout === 'undefined') ? {} : jLayout;

	jLayout.border = function(spec) {
		var my = {}, that = {}, east = spec.east, west = spec.west, north = spec.north, south = spec.south, center = spec.center;

		my.hgap = spec.hgap || 0;
		my.vgap = spec.vgap || 0;

		that.items = function() {
			var items = [];
			if (east) {
				items.push(east);
			}

			if (west) {
				items.push(west);
			}

			if (north) {
				items.push(north);
			}

			if (south) {
				items.push(south);
			}

			if (center) {
				items.push(center);
			}
			return items;
		};

		that.layout = function(container) {
			var size = container.bounds(), insets = container.insets(), top = insets.top, bottom = size.height
					- insets.bottom, left = insets.left, right = size.width
					- insets.right, tmp;

			if (north && north.isVisible()) {
				tmp = north.preferredSize();
				north.bounds({
					'x' : left,
					'y' : top,
					'width' : right - left,
					'height' : tmp.height
				});
				north.doLayout();

				top += tmp.height + my.vgap;
			}
			if (south && south.isVisible()) {
				tmp = south.preferredSize();
				south.bounds({
					'x' : left,
					'y' : bottom - tmp.height,
					'width' : right - left,
					'height' : tmp.height
				});
				south.doLayout();

				bottom -= tmp.height + my.vgap;
			}
			if (east && east.isVisible()) {
				tmp = east.preferredSize();
				east.bounds({
					'x' : right - tmp.width,
					'y' : top,
					'width' : tmp.width,
					'height' : bottom - top
				});
				east.doLayout();

				right -= tmp.width + my.hgap;
			}
			if (west && west.isVisible()) {
				tmp = west.preferredSize();
				west.bounds({
					'x' : left,
					'y' : top,
					'width' : tmp.width,
					'height' : bottom - top
				});
				west.doLayout();

				left += tmp.width + my.hgap;
			}
			if (center && center.isVisible()) {
				center.bounds({
					'x' : left,
					'y' : top,
					'width' : right - left,
					'height' : bottom - top
				});
				center.doLayout();
			}
			return container;
		};

		function typeLayout(type) {
			return function(container) {
				var insets = container.insets(), width = 0, height = 0, type_size;

				if (east && east.isVisible()) {
					type_size = east[type + 'Size']();
					width += type_size.width + my.hgap;
					height = type_size.height;
				}
				if (west && west.isVisible()) {
					type_size = west[type + 'Size']();
					width += type_size.width + my.hgap;
					height = Math.max(type_size.height, height);
				}
				if (center && center.isVisible()) {
					type_size = center[type + 'Size']();
					width += type_size.width;
					height = Math.max(type_size.height, height);
				}
				if (north && north.isVisible()) {
					type_size = north[type + 'Size']();
					width = Math.max(type_size.width, width);
					height += type_size.height + my.vgap;
				}
				if (south && south.isVisible()) {
					type_size = south[type + 'Size']();
					width = Math.max(type_size.width, width);
					height += type_size.height + my.vgap;
				}

				return {
//					'width' : width + insets.left + insets.right,
					'width' : 'auto',
					'height' : height + insets.top + insets.bottom
				};
			};
		}
		that.preferred = typeLayout('preferred');
		that.minimum = typeLayout('minimum');
		that.maximum = typeLayout('maximum');
		return that;
	};
}());

if (jQuery && jLayout) {
	(function($) {
		/**
		 * This wraps jQuery objects in another object that supplies the methods
		 * required for the layout algorithms.
		 */
		function wrap(item, resize) {
			var that = {};

			$.each([ 'min', 'max' ], function(i, name) {
				that[name + 'imumSize'] = function(value) {
					var l = item.data('jlayout');

					if (l) {
						return l[name + 'imum'](that);
					} else {
						return item[name + 'Size'](value);
					}
				};
			});

			$
					.extend(
							that,
							{
								doLayout : function() {
									var l = item.data('jlayout');

									if (l) {
										l.layout(that);
									}
									item.css({
										position : 'absolute'
									});
								},
								isVisible : function() {
									return item.isVisible();
								},
								insets : function() {
									var p = item.padding(), b = item.border();

									return {
										'top' : p.top,
										'bottom' : p.bottom + b.bottom + b.top,
										'left' : p.left,
										'right' : p.right + b.right + b.left
									};
								},
								bounds : function(value) {
									var tmp = {};

									if (value) {
										if (typeof value.x === 'number') {
											tmp.left = value.x;
										}
										if (typeof value.y === 'number') {
											tmp.top = value.y;
										}
										if (typeof value.width === 'number') {
											tmp.width = (value.width - (item
													.outerWidth(true) - item
													.width()));
											tmp.width = (tmp.width >= 0) ? tmp.width
													: 0;
										}
										if (typeof value.height === 'number') {
											tmp.height = value.height
													- (item.outerHeight(true) - item
															.height());
											tmp.height = (tmp.height >= 0) ? tmp.height
													: 0;
										}
										item.css(tmp);
										return item;
									} else {
										tmp = item.position();
										return {
											'x' : tmp.left,
											'y' : tmp.top,
											'width' : item.outerWidth(false),
											'height' : item.outerHeight(false)
										};
									}
								},
								preferredSize : function() {
									var minSize, maxSize, margin = item
											.margin(), size = {
										width : 0,
										height : 0
									}, l = item.data('jlayout');

									if (l && resize) {
										size = l.preferred(that);

										minSize = that.minimumSize();
										maxSize = that.maximumSize();

										size.width += margin.left
												+ margin.right;
										size.height += margin.top
												+ margin.bottom;

										if (size.width < minSize.width
												|| size.height < minSize.height) {
											size.width = Math.max(size.width,
													minSize.width);
											size.height = Math.max(size.height,
													minSize.height);
										} else if (size.width > maxSize.width
												|| size.height > maxSize.height) {
											size.width = Math.min(size.width,
													maxSize.width);
											size.height = Math.min(size.height,
													maxSize.height);
										}
									} else {
										size = that.bounds();
										size.width += margin.left
												+ margin.right;
										size.height += margin.top
												+ margin.bottom;
									}
									return size;
								}
							});
			return that;
		}

		$.fn.layout = function(options) {
			var opts = $.extend({}, $.fn.layout.defaults, options);
			return $
					.each(
							this,
							function() {
								var element = $(this), o = (element
										.data('layout')) ? $.extend(opts,
										element.data('layout')) : opts, elementWrapper = wrap(
										element, o.resize);

								if (o.type === 'border'
										&& typeof jLayout.border !== 'undefined') {
									$.each([ 'north', 'south', 'west', 'east',
											'center' ], function(i, name) {
										if (element.children().hasClass(name)) {
											o[name] = wrap(element.children('.'
													+ name + ':first'));
										}
									});
									element.data('jlayout', jLayout.border(o));
								} else if (o.type === 'grid'
										&& typeof jLayout.grid !== 'undefined') {
									o.items = [];
									element
											.children()
											.each(
													function(i) {
														if (!$(this)
																.hasClass(
																		'ui-resizable-handle')) {
															o.items[i] = wrap($(this));
														}
													});
									element.data('jlayout', jLayout.grid(o));
								} else if (o.type === 'flexGrid'
										&& typeof jLayout.flexGrid !== 'undefined') {
									o.items = [];
									element
											.children()
											.each(
													function(i) {
														if (!$(this)
																.hasClass(
																		'ui-resizable-handle')) {
															o.items[i] = wrap($(this));
														}
													});
									element
											.data('jlayout', jLayout
													.flexGrid(o));
								} else if (o.type === 'column'
										&& typeof jLayout.column !== 'undefined') {
									o.items = [];
									element
											.children()
											.each(
													function(i) {
														if (!$(this)
																.hasClass(
																		'ui-resizable-handle')) {
															o.items[i] = wrap($(this));
														}
													});
									element.data('jlayout', jLayout.column(o));
								} else if (o.type === 'flow'
										&& typeof jLayout.flow !== 'undefined') {
									o.items = [];
									element
											.children()
											.each(
													function(i) {
														if (!$(this)
																.hasClass(
																		'ui-resizable-handle')) {
															o.items[i] = wrap($(this));
														}
													});
									element.data('jlayout', jLayout.flow(o));
								}

								if (o.resize) {
									elementWrapper.bounds(elementWrapper
											.preferredSize());
								}

								elementWrapper.doLayout();
								element.css({
									position : 'relative'
								});
								if ($.ui !== undefined) {
									element.addClass('ui-widget');
								}
							});
		};

		$.fn.layout.defaults = {
			resize : true,
			type : 'grid'
		};
	}(jQuery));
}

jQuery(function($) {
	var container = $('.layout');

	function layout() {
		container.layout({
			resize : false,
			type : 'border',
			vgap : 8,
			hgap : 8
		});
	}

	 $('#north').resizable({
	 handles : 's',
	 stop : layout,
	 helper : 'ui-resizable-helper-north'
	 });
	 
//	 $('.north').resizable({
//	 handles : 's',
//	 stop : layout,
//	 helper : 'ui-resizable-helper-north'
//	 });

	$('.south').resizable({
		handles : 'n',
		stop : layout,
		helper : 'ui-resizable-helper-south'
	});

	$('.east').resizable({
		handles : 'w',
		stop : layout,
		helper : 'ui-resizable-helper-east'
	});

	$('.west').resizable({
		handles : 'e',
		stop : layout,
		helper : 'ui-resizable-helper-west'
	});

	$(window).resize(layout);

	layout();
});