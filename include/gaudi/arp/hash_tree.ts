// @ts-nocheck

const UNULL = -1;

/*
references: 
https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/
*/

function sgn(val) { return val < 0 ? -1 : 1; }

function dump_binary(i) {
  let s = "";
  for (let j = 31; j >= 0; j--) {
    s += ((i & (1 << j)) ? '1' : '0');
  }
  return s;
}

function scale(x, c) {
  return Math.min(Math.max(x * c, 0.0), c - 1.0);
}
const magicBitsMask2DEncode = [
  0xffffffff, 0x0000ffff, 0x00ff00ff, 0x0f0f0f0f, 0x33333333, 0x55555555,
];

function morton2DSplitBy2bits(coord) {
  const masks = magicBitsMask2DEncode;
  let x = coord & masks[0];
  x = (x | (x << 16)) & masks[1];
  x = (x | (x << 8)) & masks[2];
  x = (x | (x << 4)) & masks[3];
  x = (x | (x << 2)) & masks[4];
  x = (x | (x << 1)) & masks[5];
  return x;
}

export function morton2DEncode(x, y) {
  const max15 = Math.pow(2, 15) - 1;

  console.assert(x >= 0 && x <= 1);
  console.assert(y >= 0 && y <= 1);
  x = Math.floor(x * max15);
  y = Math.floor(y * max15);
  //if (!valuesAreUint15(x, y)) {
  //  throw CoordRange15Error;
  //}
  return morton2DSplitBy2bits(x) | (morton2DSplitBy2bits(y) << 1);
}

function clz(i, j, hash) {

  if (j < 0)
    return -1;
  if (j > hash.length - 1)
    return -1;
  const code_i = hash[i];
  const code_j = hash[j];
  /*
  console.log(i, j);
  console.log(dump_binary(code_i ^ code_j));
  console.log(decodeMorton2(code_i));
  console.log(decodeMorton2(code_j));
  */
  return Math.clz32(code_i ^ code_j);
}

function find_split(start, end, hash) {

  // uint first_code = hash[ids[start]];
  // uint last_code = hash[ids[end]];
  // if (first_code == last_code)
  //   return (first_code + last_code) >> 1;

  // uint common_prefix = std::countl_zero(first^last);

  const common_prefix_dist = clz(start, end, hash);
  let split = start;
  let step = end - start;
  while (step > 1) {
    step = (step + 1) >> 1; // exponential decrease
    let new_split = split + step;
    if (new_split < end) {

      const new_prefix_dist = clz(start, new_split, hash);

      if (new_prefix_dist > common_prefix_dist) {
        split = new_split;
      }
    }
  }

  return split;
}

function find_range(i, hash) {
  const N = hash.length;

  const dir = sgn(clz(i, i + 1, hash) - clz(i, i - 1, hash));
  const sig_min = clz(i, i - dir, hash);

  let lmax = 2;
  while (clz(i, i + lmax * dir, hash) > sig_min)
    lmax *= 2;

  let l = 0;
  let t = lmax;

  while (t >= 1) {
    t = t >> 1;
    if (clz(i, i + (l + t) * dir, hash) > sig_min)
      l += t;
  }
  const j = i + l * dir;

  return dir < 0 ? [j, i] : [i, j];
}

class radix_tree_node {
  start = UNULL;
  end = UNULL;
  split = UNULL; // split is the index of the last element in the left child
  parent = UNULL;
};

function dump_nodes(nodes) {
  nodes.forEach((element, i) => {
    console.log("node[", i, "]: ", element.start, element.end, element.split, element.parent);
  });
}


export function build_tree(ids, hash) {
  let nodes = Array.from({ length: ids.length + ids.length - 1 }, (_, i) => new radix_tree_node());
  //nodes = nodes.map(e => new radix_tree_node());

  const leaf_start = ids.length - 1;
  for (let i = 0; i < ids.length; i++) {
    nodes[leaf_start + i].start = leaf_start + i;
    nodes[leaf_start + i].end = leaf_start + i + 1;
    nodes[leaf_start + i].split = UNULL;
    nodes[leaf_start + i].parent = UNULL;
  }

  for (let i = 0; i < ids.length - 1; i++) {
    const range = find_range(i, hash);
    const split = find_split(range[0], range[1], hash);

    nodes[i].start = range[0];
    nodes[i].end = range[1];
    nodes[i].split = split;
    if (split == range[0])
      nodes[leaf_start + range[0]].parent = i;
    else
      nodes[split].parent = i;

    if (split + 1 == range[1])
      nodes[leaf_start + range[1]].parent = i;
    else
      nodes[split + 1].parent = i;
  }

  if (0) // dump nodes
    dump_nodes(nodes);

  return nodes;
}

export function test_tree(nodes, ids, hash) {
  let visited = Array.from({ length: ids.length }, (_, i) => false);

  let stack = [];
  stack.push(0);
  let leaf_startt = ids.length - 1;
  const leaf_start = (nodes.length + 1) / 2 - 1;
  let no_asserts = true;
  while (stack.length > 0) {
    const i = stack.pop();
    const lid = i - leaf_start;
    //console.log("node[", i, "]", nodes[i].start, nodes[i].end, nodes[i].split, nodes[i].parent);

    if (nodes[i].split + 0 == nodes[i].start) {
      const j0 = nodes[i].start;
      visited[lid] = true;
      console.assert(nodes[leaf_start + j0].parent == i);

    } else if (nodes[i].split + 1 == nodes[i].end) {
      const j1 = nodes[i].end;
      visited[lid] = true;
      console.assert(nodes[leaf_start + j1].parent == i);
    }
    else {
      let start_i = nodes[i].start;
      let end_i = nodes[i].end;

      let j0 = nodes[i].split;
      let j1 = nodes[i].split + 1;
      // assert parents
      console.assert(nodes[j0].parent == i);
      console.assert(nodes[j1].parent == i);
      // assert rangings
      console.assert(nodes[j0].start >= nodes[i].start);
      console.assert(nodes[j0].end <= nodes[i].end);
      console.assert(nodes[j1].start >= nodes[i].start);
      console.assert(nodes[j1].end <= nodes[i].end);
      stack.push(j0);
      stack.push(j1);
    }
  }
  let visited_all = false;
  for (let i = 0; i < visited.length; i++) {
    if (!visited[i]) {
      visited_all = false;
      break;
    }
  }
}

export function unit_test_tree() {
  const hash = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
  // std::vector<int> hash = {0, 1, 2, 3, 4, 5, 6, 7};
  const ids = hash;
  const nodes = build_tree(ids, hash);
  test_tree(nodes, ids, hash);
}

// function object prototype that takes T A and B and returns T, to generically perform pyramid ops
const is_leaf = (node) => {
  return node.split == node.start || node.split + 1 == node.end;
}

const get_leaf_idx = () => {
  //ids.length + ids.length - 1
}

export function traverse(nodes, fcn) {
  if (nodes.length === 0) return;
  let stack = [0];
  while (stack.length > 0) {
    let cid = stack.pop();
    let cnode = nodes[cid];
    if (cnode.split == UNULL) continue;

    const leaf_start = (nodes.length + 1) / 2 - 1;
    if (nodes[cid].split + 0 == nodes[cid].start &&
      nodes[cid].split + 1 == nodes[cid].end) {
      const j0 = nodes[cid].start;
      const j1 = nodes[cid].end;
      fcn(cid, j0, nodes[leaf_start + j0])
      fcn(cid, j1, nodes[leaf_start + j1])
    } else if (nodes[cid].start == nodes[cid].end) {
      /*why is this node getting generated?*/
      const j0 = nodes[cid].start;
      fcn(cid, j0, nodes[leaf_start + j0]);
    }
    else if (nodes[cid].split + 0 == nodes[cid].start) {
      const j0 = nodes[cid].start;
      fcn(cid, j0, nodes[leaf_start + j0]);
      stack.push(nodes[cid].split + 1);
    } else if (nodes[cid].split + 1 == nodes[cid].end) {
      const j1 = nodes[cid].end;
      fcn(cid, j1, nodes[leaf_start + j1]);
      stack.push(nodes[cid].split);
    } else if (fcn(cid, -1, cnode)) {
      //console.log(cid, lid, cnode);
      stack.push(nodes[cid].split);
      stack.push(nodes[cid].split + 1);
    };
  }

}

export function build_pyramid(data, idx, nodes, in_fnc, op_fnc) {
  // THE GPU is probably going to need a max tree depth, but we can just do a max recursion rate of 32?
  const max_depth = 32;
  let pyramid = Array.from({ length: nodes.length }, in_fnc);
  const leaf_start = data.length - 1;
  for (let i = 0; i < data.length; i++) {
    // insert leaf then propagate up the pyramid

    const ii = leaf_start + i;


    let data_i = data[idx[i]];
    pyramid[ii] = data_i; //set leaf nodes to data

    let parent = nodes[ii].parent;

    let j = 0;
    while (j < max_depth && parent != UNULL) {

      pyramid[parent] = op_fnc(pyramid[parent], data_i);
      data_i = pyramid[parent];
      parent = nodes[parent].parent;
      j++;
    }
  }
  return pyramid;
}

function normalize_data(data) {
  let min = [Infinity, Infinity, Infinity];
  let max = [-Infinity, -Infinity, -Infinity];
  let normalizedData = data.map(e => e.slice());

  for (let i = 0; i < data.length; i++) {
    for (let j = 0; j < 3; j++) {
      min[j] = Math.min(min[j], data[i][j]);
      max[j] = Math.max(max[j], data[i][j]);
    }
  }
  for (let i = 0; i < normalizedData.length; i++) {
    for (let j = 0; j < 3; j++) {
      let l = max[j] - min[j];
      l = l == 0 ? 1 : l;
      normalizedData[i][j] = (data[i][j] - min[j]) / l;
    }
  }
  return normalizedData;
}

export function mk_hash(data) {
  const normalizedData = normalize_data(data);

  const hashes = normalizedData.map((e, i) => morton2DEncode(e[0], e[2]));
  //const hashes = normalizedData.map((e, i) => EncodeMorton3(e[0], e[1], e[2]));
  const indexes = Array.from({ length: hashes.length }, (_, i) => i);
  const sortedIndices = indexes.sort((a, b) => hashes[a] - hashes[b]);
  const sortedHashes = sortedIndices.map(i => hashes[i]);
  return [sortedHashes, sortedIndices];
}

