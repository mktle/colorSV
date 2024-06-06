#!/usr/bin/env k8

/**************
 * Author: Heng Li *
 **************/

const gc_version = "r115";

/**************
 * From k8.js *
 **************/

Array.prototype.delete_at = function(i) {
	for (let j = i; j < this.length - 1; ++j)
		this[j] = this[j + 1];
	--this.length;
}

function* getopt(argv, ostr, longopts) {
	if (argv.length == 0) return;
	let pos = 0, cur = 0;
	while (cur < argv.length) {
		let lopt = "", opt = "?", arg = "";
		while (cur < argv.length) { // skip non-option arguments
			if (argv[cur][0] == "-" && argv[cur].length > 1) {
				if (argv[cur] == "--") cur = argv.length;
				break;
			} else ++cur;
		}
		if (cur == argv.length) break;
		let a = argv[cur];
		if (a[0] == "-" && a[1] == "-") { // a long option
			pos = -1;
			let c = 0, k = -1, tmp = "", o;
			const pos_eq = a.indexOf("=");
			if (pos_eq > 0) {
				o = a.substring(2, pos_eq);
				arg = a.substring(pos_eq + 1);
			} else o = a.substring(2);
			for (let i = 0; i < longopts.length; ++i) {
				let y = longopts[i];
				if (y[y.length - 1] == "=") y = y.substring(0, y.length - 1);
				if (o.length <= y.length && o == y.substring(0, o.length)) {
					k = i, tmp = y;
					++c; // c is the number of matches
					if (o == y) { // exact match
						c = 1;
						break;
					}
				}
			}
			if (c == 1) { // find a unique match
				lopt = tmp;
				if (pos_eq < 0 && longopts[k][longopts[k].length-1] == "=" && cur + 1 < argv.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				}
			}
		} else { // a short option
			if (pos == 0) pos = 1;
			opt = a[pos++];
			let k = ostr.indexOf(opt);
			if (k < 0) {
				opt = "?";
			} else if (k + 1 < ostr.length && ostr[k+1] == ":") { // requiring an argument
				if (pos >= a.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				} else arg = a.substring(pos);
				pos = -1;
			}
		}
		if (pos < 0 || pos >= argv[cur].length) {
			argv.delete_at(cur);
			pos = 0;
		}
		if (lopt != "") yield { opt: `--${lopt}`, arg: arg };
		else if (opt != "?") yield { opt: `-${opt}`, arg: arg };
		else yield { opt: "?", arg: "" };
	}
}

function* k8_readline(fn) {
	let buf = new Bytes();
	let file = new File(fn);
	while (file.readline(buf) >= 0) {
		yield buf.toString();
	}
	file.close();
	buf.destroy();
}

function parseNum(s) {
	var m, x = null;
	if ((m = /^(\d*\.?\d*)([mMgGkK]?)/.exec(s)) != null) {
		x = parseFloat(m[1]);
		if (m[2] == 'k' || m[2] == 'K') x *= 1000;
		else if (m[2] == 'm' || m[2] == 'M') x *= 1000000;
		else if (m[2] == 'g' || m[2] == 'G') x *= 1000000000;
	}
	return Math.floor(x + .499);
}

/********************************
 * Extract SVs from GAF/PAF/SAM *
 ********************************/

function mg_revcomp(s) {
	function complement(x) { return { a:'t', t:'a', g:'c', c:'g' }[x] }
	return s.split('').reverse().map(complement).join('');
}

function gc_cmd_extract(args) {
	let opt = { min_mapq:5, min_mapq_end:30, min_frac:0.7, min_len:100, min_aln_len_end:2000, min_aln_len_mid:50, max_cnt_10k:3,
				dbg:false, polyA_pen:5, polyA_drop:100, name:"foo", cen:{} };
	for (const o of getopt(args, "q:Q:l:dc:a:e:m:n:b:", [])) {
		if (o.opt == "-q") opt.min_mapq = parseInt(o.arg);
		else if (o.opt == "-Q") opt.min_mapq_end = parseInt(o.arg);
		else if (o.opt == "-l") opt.min_len = parseNum(o.arg);
		else if (o.opt == "-d") opt.dbg = true;
		else if (o.opt == "-f") opt.min_frac = parseFloat(o.arg);
		else if (o.opt == "-c") opt.max_cnt_10k = parseInt(o.arg);
		else if (o.opt == "-a") opt.polyA_pen = parseInt(o.arg);
		else if (o.opt == "-e") opt.min_aln_len_end = parseInt(o.arg);
		else if (o.opt == "-m") opt.min_aln_len_mid = parseInt(o.arg);
		else if (o.opt == "-n") opt.name = o.arg;
		else if (o.opt == "-b") {
			for (const line of k8_readline(o.arg)) {
				const t = line.split("\t");
				if (opt.cen[t[0]] == null) opt.cen[t[0]] = [];
				opt.cen[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
			}
			for (const ctg in opt.cen)
				opt.cen[ctg].sort(function(a,b) { return a[0]-b[0] });
		}
	}
	if (opt.min_mapq > opt.min_mapq_end)
		opt.min_mapq = opt.min_mapq_end;
	if (args.length == 0) {
		print("Usage: gafcall.js extract [options] <stable.gaf>");
		print("Options:");
		print(`  -n STR     sample name [${opt.name}]`);
		print(`  -l INT     min INDEL len [${opt.min_len}]`);
		print(`  -f FLOAT   min mapped query fraction [${opt.min_frac}]`);
		print(`  -c INT     max number of long INDELs per 10kb [${opt.max_cnt_10k}]`);
		print(`  -q INT     min mapq [${opt.min_mapq}]`);
		print(`  -Q INT     min mapq for alignment ends [${opt.min_mapq_end}]`);
		print(`  -e INT     min alignment length at ends [${opt.min_aln_len_end}]`);
		print(`  -m INT     min alignment length in the middle [${opt.min_aln_len_mid}]`);
		print(`  -a INT     penalty for non-polyA bases [${opt.polyA_pen}]`);
		print(`  -b FILE    BED for centromeres []`);
		return;
	}

	let re = /(\d+)([=XIDMSHN])/g; // regex for CIGAR
	let re_path = /([><])([^><:\s]+):(\d+)-(\d+)/g; // regex for path/ctg
	let re_ds = /([\+\-\*:])([A-Za-z\[\]0-9]+)/g; // regex for the ds tag
	let re_tsd = /(\[([A-Za-z]+)\])?([A-Za-z]+)(\[([A-Za-z]+)\])?/; // regex for parsing TSD
	let global_qname = "N/A";

	function cal_cen_dist(opt, ctg, pos) {
		if (opt.cen[ctg] == null) return 1e9;
		let min = 1e9;
		for (let i = 0; i < opt.cen[ctg].length; ++i) { // TODO: binary search would be better
			const b = opt.cen[ctg][i];
			const d = pos < b[0]? b[0] - pos : pos < b[1]? 0 : pos - b[1];
			min = min < d? min : d;
		}
		return min;
	}

	function cal_cen_overlap(opt, ctg, st0, en0) {
		if (opt.cen[ctg] == null) return 0;
		let cov_st = 0, cov_en = 0, cov = 0;
		for (let i = 0; i < opt.cen[ctg].length; ++i) { // TODO: binary search would be better
			const b = opt.cen[ctg][i];
			if (b[1] <= st0 || b[0] >= en0) continue; // not overlapping with [st0, en0)
			const st1 = b[0] > st0? b[0] : st0;
			const en1 = b[1] < en0? b[1] : en0;
			if (st1 > cov_en) {
				cov += cov_en - cov_st;
				cov_st = st1, cov_en = en1;
			} else cov_en = cov_en > en1? cov_en : en1;
		}
		cov += cov_en - cov_st;
		return cov;
	}

	/*******************************************
	 * Extract long indels contained in CIGARs *
	 *******************************************/

	function cal_polyA_len(opt, int_seq) {
		let polyA_len = 0, polyT_len = 0, polyA_max = 0, polyT_max = 0;
		let score, max, max_j;
		// look for polyA on the 3'-end
		score = max = 0, max_j = int_seq.length;
		for (let j = int_seq.length - 1; j >= 0; --j) {
			if (int_seq[j] == 'A' || int_seq[j] == 'a') ++score;
			else score -= opt.polyA_pen;
			if (score > max) max = score, max_j = j;
			else if (max - score > opt.polyA_drop) break;
		}
		polyA_len = int_seq.length - max_j, polyA_max = max;
		// look for polyT on the 5'-end
		score = max = 0, max_j = -1;
		for (let j = 0; j < int_seq.length; ++j) {
			if (int_seq[j] == 'T' || int_seq[j] == 't') ++score;
			else score -= opt.polyA_pen;
			if (score > max) max = score, max_j = j;
			else if (max - score > opt.polyA_drop) break;
		}
		polyT_len = max_j + 1, polyT_max = max;
		// choose the longer one
		return polyA_max >= polyT_max? polyA_len : -polyT_len;
	}

	function path2ctg(seg, path_off, is_end) {
		let b = [];
		for (let i = 0, k = 0; i < path_off.length; ++i) {
			if (is_end) {
				while (k < seg.length && seg[k].path_en < path_off[i]) ++k;
			} else {
				while (k < seg.length && seg[k].path_en <= path_off[i]) ++k;
			}
			if (k == seg.length) throw Error(`failed to convert path offset to contig offset for read ${global_qname}`);
			const l = path_off[i] - seg[k].path_st;
			if (seg[k].strand > 0)
				b.push({ seg:k, pos:seg[k].ctg_st + l });
			else
				b.push({ seg:k, pos:seg[k].ctg_en - l });
		}
		return b;
	}

	function get_indel(opt, z) {
		if (z.length == 0) return;
		for (let j = 0; j < z.length; ++j) {
			const y = z[j];
			if (y.qen - y.qst < y.qlen * opt.min_frac) continue; // ignore short alignments
			const is_rev = (y.strand === "-");
			let m, a = [], x = y.tst, q = 0;
			while ((m = re.exec(y.cg)) != null) { // collect the list of long indels
				const op = m[2], len = parseInt(m[1]);
				if (len >= opt.min_len) {
					if (op === "I") {
						const qoff = is_rev? y.qen - (q + len) : q + y.qst;
						a.push({ st:x, en:x,     len:len,  indel_seq:".", tsd_len:0, tsd_seq:".", polyA_len:0, int_seq:".", qoff:qoff, qoff_l:qoff, qoff_r:qoff+len });
					} else if (op === "D") {
						const qoff = is_rev? y.qen - q : q + y.qst;
						a.push({ st:x, en:x+len, len:-len, indel_seq:".", tsd_len:0, tsd_seq:".", polyA_len:0, int_seq:".", qoff:qoff, qoff_l:qoff, qoff_r:qoff });
					}
				}
				if (op == "M" || op == "=" || op == "X" || op == "D" || op === "N")
					x += len;
				if (op == "M" || op == "=" || op == "X" || op == "I" || op === "S" || op === "H")
					q += len;
			}
			if (a.length == 0 || a.length > y.qlen * 1e-4 * opt.max_cnt_10k) continue;
			// set stl/enl and str/enr
			for (let i = 0; i < a.length; ++i) {
				a[i].stl = a[i].str = a[i].st;
				a[i].enl = a[i].enr = a[i].en;
			}
			// parse ds:Z
			if (y.ds) { // this MUST match CIGAR parsing
				let i = 0, x = y.tst, m;
				while ((m = re_ds.exec(y.ds)) != null) {
					const op = m[1], str = m[2];
					const seq = op === "+" || op === "-"? str.replace(/[\[\]]/g, "") : "";
					const len = op === ":"? parseInt(str) : op === "*"? 1 : op === "+" || op === "-"? seq.length : -1;
					if (len < 0) throw Error("can't determine length from the ds tag");
					if (len >= opt.min_len) { // extract INDEL sequence and check consistency with CIGAR
						if (op === "+") {
							if (a[i].st != x || a[i].en != x || a[i].len != len)
								throw Error(`inconsistent insertion at line ${lineno}`);
							a[i++].indel_seq = str;
						} else if (op === "-") {
							if (a[i].st != x || a[i].en != x + len || a[i].len != -len)
								throw Error(`inconsistent deletion at line ${lineno}`);
							a[i++].indel_seq = str;
						}
					}
					if (op === "*" || op === ":" || op === "-")
						x += len;
				}
				for (let i = 0; i < a.length; ++i) { // compute TSD and polyA lengths
					if ((m = re_tsd.exec(a[i].indel_seq)) == null)
						throw Error("Bug!");
					const tsd = (m[5]? m[5] : "") + (m[2]? m[2] : "");
					const int_seq = m[3]; // internal sequence
					a[i].tsd_len = tsd.length;
					a[i].tsd_seq = tsd;
					a[i].int_seq = int_seq;
					if (int_seq.length > 0)
						a[i].polyA_len = cal_polyA_len(opt, int_seq);
					const llen = m[2]? m[2].length : 0;
					const rlen = m[5]? m[5].length : 0;
					a[i].stl = a[i].st - rlen;
					a[i].enl = a[i].en - rlen;
					a[i].str = a[i].st + llen;
					a[i].enr = a[i].en + llen;
					if (is_rev) {
						a[i].qoff_l -= llen;
						a[i].qoff_r += rlen;
					} else {
						a[i].qoff_l -= rlen;
						a[i].qoff_r += llen;
					}
				}
			} // ~if(y.ds)
			if (opt.dbg) print('X0', line);
			let seg = []; // reference segments in the path
			if (/[><]/.test(y.path)) { // with ><: this is a path
				let x = 0;
				if (y.strand != '+') throw Error("reverse strand on path");
				while ((m = re_path.exec(y.path)) != null) {
					const st = parseInt(m[3]), en = parseInt(m[4]);
					const strand = m[1] === '>'? 1 : -1;
					seg.push({ ctg:m[2], ctg_st:st, ctg_en:en, strand:strand, path_st:x, path_en:x + (en - st) });
					x += en - st;
				}
			} else { // this is a contig name
				seg.push({ ctg:y.path, ctg_st:0, ctg_en:y.tlen, strand:1, path_st:0, path_en:y.tlen });
			}
			let off_stl = [], off_str = [], off_enl = [], off_enr = [];
			for (let i = 0; i < a.length; ++i) {
				off_stl[i] = a[i].stl, off_enl[i] = a[i].enl;
				off_str[i] = a[i].str, off_enr[i] = a[i].enr;
			}
			global_qname = y.qname;
			const stl = path2ctg(seg, off_stl, false);
			const enl = path2ctg(seg, off_enl, true);
			const str = path2ctg(seg, off_str, false);
			const enr = path2ctg(seg, off_enr, true);
			for (let i = 0; i < a.length; ++i) {
				if (!(stl[i].seg == str[i].seg && stl[i].seg == enl[i].seg && str[i].seg == enr[i].seg)) continue; // all on the same segment
				const s = seg[stl[i].seg];
				let st = stl[i].pos, en = enl[i].pos, strand = y.strand;
				if (s.strand < 0) { // then reverse complement tsd, polyA and insert
					a[i].polyA_len = -a[i].polyA_len;
					a[i].tsd_seq = mg_revcomp(a[i].tsd_seq);
					a[i].int_seq = mg_revcomp(a[i].int_seq);
					st = enr[i].pos, en = str[i].pos;
					strand = y.strand === "+"? "-" : "+";
				}
				let info1 = (a[i].len > 0? "SVTYPE=INS" : "SVTYPE=DEL") + `;SVLEN=${a[i].len};qoff_l=${a[i].qoff_l};qoff_r=${a[i].qoff_r};tsd_len=${a[i].tsd_len};polyA_len=${a[i].polyA_len}`;
				const info2 = `source=${opt.name};tsd_seq=${a[i].tsd_seq.length>0?a[i].tsd_seq:"."};insert=${a[i].int_seq.length>0?a[i].int_seq:"."}`;
				if (opt.cen[s.ctg] != null) {
					const dist_st = cal_cen_dist(opt, s.ctg, st);
					const dist_en = cal_cen_dist(opt, s.ctg, en);
					info1 += `;cen_dist=${dist_st < dist_en? dist_st : dist_en}`
				}
				print(s.ctg, st, en, y.qname, y.mapq, strand, `${info1};${info2}`);
			} // ~for(i)
		} // ~for(j)
	} // ~get_indel()

	/*********************************
	 * Extract alignment breakpoints *
	 *********************************/

	function get_end_coor(y) {
		let r = [{}, {}];
		if (/^[><]/.test(y.path)) { // a path
			if (y.strand != '+') throw Error("reverse strand on path");
			let x = 0, i = 0;
			while ((m = re_path.exec(y.path)) != null) {
				const st = parseInt(m[3]), en = parseInt(m[4]), len = en - st;
				if (y.tst >= x && y.tst < x + len) {
					r[0] = { ctg:m[2], ori:m[1], pos:-1 };
					r[0].pos = m[1] === ">"? st + (y.tst - x) : st + (x + len - y.tst) - 1;
				}
				if (y.ten > x && y.ten <= x + len) {
					r[1] = { ctg:m[2], ori:m[1], pos:-1 };
					r[1].pos = m[1] === ">"? st + (y.ten - x) - 1 : st + (x + len - y.ten);
				}
				x += len;
			}
		} else { // a contig
			r[0] = { ctg:y.path, ori: y.strand === "+"? ">" : "<", pos:-1 };
			r[0].pos = y.strand === "+"? y.tst : y.ten - 1;
			r[1] = { ctg:y.path, ori: y.strand === "+"? ">" : "<", pos:-1 };
			r[1].pos = y.strand === "+"? y.ten - 1 : y.tst;
		}
		r[0].ql = r[1].ql = y.qen - y.qst;
		return r;
	}

	function infer_svtype(opt, c0, c1, ori, qgap) { // NB: c0 MUST have the smaller coordinate
		if (c0.ctg != c1.ctg) return { st:-1, en:-1, str:"SVTYPE=BND" };
		const l = c1.pos - c0.pos + 1;
		if (l < 0) throw Error("Bug!");
		if (ori === ">>" && qgap < l && l - qgap >= opt.min_len) { // deletion
			const st = qgap < 0? c0.pos + qgap : c0.pos;
			const en = qgap < 0? c1.pos + 1 - qgap : c1.pos + 1;
			return { st:st, en:en, str:`SVTYPE=DEL;SVLEN=${-(l - qgap)};sv_region=${st},${en};tsd_len=${qgap < 0? -qgap : 0}` };
		}
		if (ori === ">>" && l < qgap && qgap - l >= opt.min_len) // insertion without TSD
			return { st:c0.pos, en:c1.pos+1, str:`SVTYPE=INS;SVLEN=${qgap - l};sv_region=${c0.pos},${c1.pos+1}` };
		if (ori === "<<" && qgap > 0 && (l < c0.ql || l < c1.ql) && qgap + l >= opt.min_len) // insertion with TSD
			return { st:c0.pos, en:c1.pos+1, str:`SVTYPE=INS;SVLEN=${qgap + l};sv_region=${c0.pos},${c1.pos+1};tsd_len=${l}` }; // TODO: is sv_region correct?
		if (ori === "<<" && qgap + l >= opt.min_len) { // tandem duplication; similar to insertion with TSD
			const st = qgap < 0? c0.pos : c0.pos > qgap? c0.pos - qgap : 0;
			const en = qgap < 0? c1.pos + 1 : c1.pos + 1 + qgap;
			return { st:st, en:en, str:`SVTYPE=DUP;SVLEN=${qgap + l};sv_region=${st},${en}` };
		}
		if ((ori === "<>" || ori === "><") && l >= opt.min_len) { // inversion
			const st = qgap < 0? c0.pos + qgap : c0.pos;
			const en = qgap < 0? c1.pos + 1 - qgap : c1.pos + 1;
			return { st:st, en:en, str:`SVTYPE=INV;SVLEN=${l - qgap};sv_region=${st},${en}` };
		}
		return { st:-1, en:-1, str:"SVTYPE=BND" };
	}

	function get_breakpoint(opt, z) {
		if (z.length < 2) return;
		z.sort(function(a,b) { return a.qst - b.qst }); // sort by start position on the read
		// filter out short alignment towards the end of the read
		let zen = z.length;
		for (let j = z.length - 1; j >= 0; --j) {
			const y = z[j];
			if (y.qen - y.qst < opt.min_aln_len_end || y.mapq < opt.min_mapq_end) zen = j;
			else break;
		}
		if (zen < 2) return;
		// filter out short alignment towards the start of the read
		let zst = 0;
		for (let j = 0; j < zen; ++j) {
			const y = z[j];
			if (y.qen - y.qst < opt.min_aln_len_end || y.mapq < opt.min_mapq_end) zst = j + 1;
			else break;
		}
		if (zen - zst < 2) return;
		// construct the final alignment list
		let zz = [];
		for (let j = zst; j < zen; ++j)
			if (z[j].qen - z[j].qst >= opt.min_aln_len_mid)
				zz.push(z[j]);
		if (zz.length < 2) return; // shouldn't happen if mid<end; just in case
		// calculate the end coordinates
		for (let j = 0; j < zz.length; ++j) {
			const r = get_end_coor(zz[j]);
			zz[j].coor = r;
		}
		// extract alignment breakpoints
		for (let j = 1; j < zz.length; ++j) {
			const y0 = zz[j-1], y1 = zz[j];
			const qgap = y1.qst - y0.qen;
			let c0 = y0.coor[1], c1 = y1.coor[0], strand2 = "+", ori = c0.ori + c1.ori;
			if (!(c0.ctg < c1.ctg || (c0.ctg === c1.ctg && c0.pos < c1.pos)))
				c0 = y1.coor[0], c1 = y0.coor[1], strand2 = "-", ori = (c1.ori === ">"? "<" : ">") + (c0.ori === ">"? "<" : ">");
			const sv_info = infer_svtype(opt, c0, c1, ori, qgap);
			let cen_str = "";
			if (opt.cen[c0.ctg] != null || opt.cen[c1.ctg] != null) {
				const dist0 = cal_cen_dist(opt, c0.ctg, c0.pos);
				const dist1 = cal_cen_dist(opt, c1.ctg, c1.pos);
				cen_str = `;cen_dist=${dist0<dist1?dist0:dist1}`;
				if (sv_info.st >= 0 && sv_info.en >= sv_info.st) {
					const ov = cal_cen_overlap(opt, c0.ctg, sv_info.st, sv_info.en);
					cen_str += `;cen_overlap=${ov}`;
				}
			}
			const qoff_l = y0.qen < y1.qst? y0.qen : y1.qst;
			const qoff_r = y0.qen > y1.qst? y0.qen : y1.qst;
			print(c0.ctg, c0.pos, ori, c1.ctg, c1.pos, y0.qname, y0.mapq < y1.mapq? y0.mapq : y1.mapq, strand2,
				  `${sv_info.str};qoff_l=${qoff_l};qoff_r=${qoff_r};qgap=${qgap};mapq=${y0.mapq},${y1.mapq};aln_len=${y0.qen-y0.qst},${y1.qen-y1.qst}${cen_str};source=${opt.name}`);
		}
	} // ~get_breakpoint()

	let lineno = 0, z = [];
	for (const line of k8_readline(args[0])) {
		++lineno;
		let t = line.split("\t");
		if (t.length < 11) continue; // SAM has 11 columns at least; PAF has 12 columns at least
		if (z.length > 0 && t[0] != z[0].qname) {
			get_indel(opt, z);
			get_breakpoint(opt, z);
			z = [];
		}
		// parse format
		let y = { qname:t[0], mapq:0, qst:-1, qen:-1, qlen:-1, tlen:-1, tst:-1, cg:null, ds:null, path:null, strand:null };
		if (t.length >= 12 && (t[4] === "+" || t[4] === "-")) { // parse PAF or GAF
			y.mapq = parseInt(t[11]);
			if (y.mapq < opt.min_mapq) continue;
			y.qlen = parseInt(t[1]);
			y.qst = parseInt(t[2]);
			y.qen = parseInt(t[3]);
			y.strand = t[4];
			y.path = t[5];
			y.tlen = parseInt(t[6]);
			y.tst = parseInt(t[7]);
			y.ten = parseInt(t[8]);
			let tp = "";
			for (let i = 12; i < t.length; ++i) {
				if (t[i].substr(0, 5) === "cg:Z:")
					y.cg = t[i].substr(5);
				else if (t[i].substr(0, 5) === "ds:Z:")
					y.ds = t[i].substr(5);
				else if (t[i].substr(0, 5) === "tp:A:")
					tp = t[i].substr(5);
			}
			if (tp != "P") continue; // filter out secondary alignments
			if (y.cg == null) continue;
		} else { // parse SAM
			if (t[0][0] === "@") continue;
			const flag = parseInt(t[1]);
			if (flag & 0x100) continue;
			y.mapq = parseInt(t[4]);
			if (y.mapq < opt.min_mapq) continue;
			y.strand = (flag&0x10)? "-" : "+";
			y.path = t[2];
			y.tlen = 0xffffffff; // tlen doesn't need to be accurate for SAM or PAF
			y.tst = parseInt(t[3]) - 1;
			y.ten = -1;
			y.cg = t[5];
			let m;
			y.qst = (m = /^(\d+)[SH]/.exec(cg)) != null? parseInt(m[1]) : 0;
			y.qlen = 0;
			while ((m = re.exec(cg)) != null) {
				const op = m[2];
				if (op == "S" || op == "H" || op == "M" || op == "=" || op == "I")
					y.qlen += parseInt(m[1]);
			}
			for (let i = 11; i < t.length; ++i)
				if (t[i].substr(0, 5) === "ds:Z:")
					y.ds = t[i].substr(5);
			y.qen = y.qlen - ((m = /(\d+)[SH]$/.exec(cg)) != null? parseInt(m[1]) : 0);
		}
		z.push(y);
	}
	get_indel(opt, z);
	get_breakpoint(opt, z);
}

/*****************
 * Main function *
 *****************/

function main(args)
{
	var cmd = args.shift();
	if (cmd === "extract" || cmd === "getsv") gc_cmd_extract(args);
	else if (cmd === "version") {
		print(gc_version);
		return;
	} else throw Error("unrecognized command: " + cmd);
}

main(arguments);
