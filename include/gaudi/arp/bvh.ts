// @ts-nocheck


import * as THREE from 'three';
import { mk_hash, build_tree, build_pyramid, traverse, test_tree } from './hash_tree.js';
import { logger } from './logger.js';

/*
let testId = 0;
function incrementTestId(){
setTimeout(() => {
  testId+=1;
  incrementTestId();
}, 1000);
}
incrementTestId();
*/
export class bvh{
  hash = [];
  idx = [];
  tree = [];
  node_bbs = [];
  centers = [];
  constructor(bounding_boxes){
    const centers = bounding_boxes.map(bb => bb.getCenter(new THREE.Vector3()));
    const centers_arr = centers.map(v => [v.x, v.y, v.z]);
    
    [this.hash, this.idx] = mk_hash(centers_arr);
    this.tree = build_tree(this.idx, this.hash);
    this.centers = centers;
    this.node_bbs = build_pyramid(bounding_boxes, this.idx, this.tree, 
      (e,i)=> new THREE.Box3(),
      (bba, bbb)=> bba.clone().union(bbb)
    );
  }

  draw = () => {
      for(let i = 0; i < this.idx.length; i++){
      let leafStart = this.idx.length - 1;
      let stack = [leafStart + i];
      //if (i != testId) continue;
      while(stack.length > 0){
        const cid = stack.pop();
        const cnode = this.tree[cid];
        if(cnode.parent == -1) break;

        const bbcen0 = this.node_bbs[cid].getCenter(new THREE.Vector3());
        const bbcen1 = this.node_bbs[cnode.parent].getCenter(new THREE.Vector3());
        //logger.add_line(bbcen0, bbcen1, 'blue');
        logger.add_bbox(this.node_bbs[cid].min, this.node_bbs[cid].max, 'turquoise');
        stack.push(cnode.parent);
      }
    }
  }

  run_test = _ =>{
    /*
    for(let i = 1; i < this.idx.length; i++){
      const bbcen0 = this.centers[this.idx[i]].clone();
      const bbcen1 = this.centers[this.idx[i-1]].clone();
      //const bbcen0 = this.centers[i].clone();
      //const bbcen1 = this.centers[i-1].clone();
      
      logger.add_line(bbcen0, bbcen1, 'red');
    }
  
     
    for(let i = 0; i < this.idx.length; i++){
      let leafStart = this.idx.length - 1;
      let stack = [leafStart + i];
      //if (i != testId) continue;
      while(stack.length > 0){
        const cid = stack.pop();
        const cnode = this.tree[cid];
        if(cnode.parent == -1) break;

        const bbcen0 = this.node_bbs[cid].getCenter(new THREE.Vector3());
        const bbcen1 = this.node_bbs[cnode.parent].getCenter(new THREE.Vector3());
        logger.add_line(bbcen0, bbcen1, 'blue');
        //logger.add_bbox(this.node_bbs[cid].min, this.node_bbs[cid].max, 'turquoise');
        stack.push(cnode.parent);
      }
    }
        
    //this.idx.forEach((e,i) => {
      //let bbcen = this.node_bbs[i].getCenter(new THREE.Vector3());
    //});
    test_tree(this.tree, this.idx, this.hash);
    */
  } 

  find_overlapping = bb => {
    let result = [];
    //logger.add_bbox(bb.min, bb.max, 'red');
    const this_center = bb.getCenter(new THREE.Vector3());

    traverse(this.tree, (cid, lid, node) => {
      const overlap = this.node_bbs[cid].intersectsBox(bb);
      //if(overlap) debugger;

      if(overlap && lid > -1){
        const dist = this_center.distanceTo(this.centers[this.idx[lid]]);
        result.push([this.idx[lid], dist]);
      }
      /*
      if(overlap){
        logger.add_scaled_bbox(this.node_bbs[cid].min, this.node_bbs[cid].max, 0.25, 'cyan');
      }
      
      if(overlap && lid > -1){
        const bbox = this.node_bbs[cid];
        const center0 = new THREE.Vector3();
        bb.getCenter(center0);
        const center1 = new THREE.Vector3();
        bbox.getCenter(center1);
        logger.add_line(center0, center1, 'magenta');
        logger.add_scaled_bbox(bbox.min, bbox.max, 0.25, 'green');
      }
      */
      return overlap;
    });
    const sorted = result.sort((a,b) => a[1] - b[1]);
    /*
    sorted.forEach((e,i) => {
      if(i < 64)
        logger.add_line(this_center, this.centers[e[0]], 'yellow');
      else
        logger.add_line(this_center, this.centers[e[0]], 'green');
    });
    */
    return sorted.map(e => e[0]);
  }
}