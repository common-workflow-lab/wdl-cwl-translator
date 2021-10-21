'use-strict';
// https://github.com/anonyco/Javascript-Fast-Light-Map-WeakMap-Set-And-WeakSet-JS-Polyfill/blob/07595a2acb8df921c6481024218200c441f87621/mapPolyfill.src.js
// Creative Commons Zero v1.0 Universal

var keycur,
	i, len,
	k, v,
	iterable,
Mapproto = {
	//length: 0,
	'delete': function( key ){
		keycur = NaNsearch( this.k, key ); // k is for keys
		if (!~keycur) return false;
		this.k.splice(keycur, 1);
		this.v.splice(keycur, 1);
		--this.size;
		return true;
	},
	'get': function( key ){
		return this.v[NaNsearch(this.k, key)]; // automagicly returns undefined if it doesn't exist
	},
	'set': function( key, value ){
		keycur = NaNsearch(this.k, key);
		if (!~keycur) // if (keycur === -1)
			this.k[keycur = this.size++] = key;
		this.v[keycur] = value;
		return this;
	},
	'has': function( key ){
		return NaNsearch(this.k, key) > -1;
	},
	'clear': function(){
		this.k.length = this.v.length = this.size = 0;
		//return undefined
	},
	'forEach': function( Func, thisArg ){
		if(thisArg)
			Func = Func.bind(thisArg);
		var i = -1, len = this.size;
		while (++i !== len) Func(this.v[i], this.k[i], this);
	},
	'entries': function( ){
		var nextIndex = 0, that = this;
		return {
			next: function() {
				return nextIndex !== that.size ? {value: [that.k[nextIndex++], that.v[nextIndex]], done: false} : {done: true};
			}
		};
	},
	'keys': function( ){
		var nextIndex = 0, that = this;
		return {
			next: function() {
				return nextIndex !== that.size ? {value: that.k[nextIndex++], done: false} : {done: true};
			}
		};
	},
	'values': function( ){
		var nextIndex = 0, that = this;
		return {
			next: function() {
				return nextIndex !== that.size ? {value: that.v[nextIndex++], done: false} : {done: true};
			}
		};
	},
	toString: function(){return '[object Map]'}
};
function NaNsearch(arr, val){ // search that compensates for NaN indexs
	if (val === val) // if val is not NaN
		return arr.indexOf(val);

	i = 0, len = arr.length;
	// Check for the first index that is not itself (i.e. NaN)
	while (arr[i] === arr[i] && ++i !== len)
		; // do nothing
	return i;
};

// Map & WeakMap polyfill
var Map =  function(raw){
	k = this.k = [];
	v = this.v = [];
	len = 0;
	if (raw !== undefined  && raw !== null){
		iterable = Object(raw);
		// split up the data into two useable streams: one for keys (k), and one for values (v)
		i = +iterable.length;
		if (i != i) // if i is NaN
			throw new TypeError('(' + (raw.toString || iterable.toString)() + ') is not iterable');

		while (i--)
			if (iterable[i] instanceof Object){
				if (!~NaNsearch(k, iterable[i][0])) // if current is not already in the array
					k[len] = iterable[i][0], v[len++] = iterable[i][1]; // len++ increments len, but returns value before increment
			} else
				throw new TypeError('Iterator value ' + iterable[i] + ' is not an entry object');
		k.reverse();
		v.reverse();
	}
	this.size = len;
};
Map.prototype = Mapproto;
