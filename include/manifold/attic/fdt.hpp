#ifndef __FDTBIN__
#define __FDTBIN__

#include <bitset>
#include <functional>
#include <stack>

#include "conj_grad.hpp"
#include "m2Includes.h"
#include "quartic/cubic.hpp"
#include "tribox_test.hpp"

#include "TIMER.h"
#include <cmath>

namespace asawa {

template <typename SPACE, typename BLOCK_T> struct LeafBlock {
  M2_TYPEDEFS;

  LeafBlock() {}
  LeafBlock(coordinate_type cen) {}
  coordinate_type mCen;
  BLOCK_T data;
};

template <typename SPACE, int N, typename BLOCK_T> struct Block {
  M2_TYPEDEFS;

  Block(T length, coordinate_type cen)
  //:_xRes(N), _yRes(N), _zRes(N)
  {

    mRes = N;
    mxRes = N;
    myRes = N;
    mzRes = N;
    msize = N * N * N;
    mdx = length / N;
    mcenter = cen;

    mlengths = coordinate_type(length, length, length, 0.0);
  }

  Block()
  //:_xRes(N), _yRes(N), _zRes(N)
  {

    mRes = N;
    mxRes = N;
    myRes = N;
    mzRes = N;
    msize = N * N * N;
  }

  Block(Block &other) { copy(other); }

  ~Block() {}

  Block &operator=(Block &rhs) {
    copy(rhs);
    return *this;
  }

  Block &copy(Block &rhs) {

    if (this != &rhs) {
      for (int i = 0; i < N * N * N; i++) {
        mchildMask[i] = rhs.mchildMask[i];
      }
      mxRes = rhs.xRes();
      myRes = rhs.yRes();
      mzRes = rhs.zRes();
      mdx = rhs.dx();
      mcenter = rhs.center();
      mlengths = rhs.lengths();
    }
  }

  void setCenter(coordinate_type cen) { mcenter = cen; }

  void setLength(T length) {
    mdx = length / mRes;
    mlengths = coordinate_type(length, length, length, 0.0);
  }

  void clear() {}

  int xRes() { return mxRes; }
  int yRes() { return myRes; }
  int zRes() { return mzRes; }
  int size() { return mnumBlocks; }

  T dx() { return mdx; }

  coordinate_type center() { return mcenter; }

  coordinate_type center(int x, int y, int z) {
    coordinate_type cen = mcenter - 0.5 * mlengths +
                          coordinate_type(x * mdx, y * mdx, z * mdx, 0) +
                          0.5 * coordinate_type(mdx, mdx, mdx, 0);
    return cen;
  }

  coordinate_type center(int index) {
    int z = (int)(index / mRes / mRes);
    int y = (int)(index - z * mRes * mRes) / mRes;
    int x = (int)(index - y * mRes + z * mRes * mRes);
    coordinate_type cen = mcenter - 0.5 * mlengths +
                          coordinate_type(x * mdx, y * mdx, z * mdx, 0) +
                          0.5 * coordinate_type(mdx, mdx, mdx, 0);
    return cen;
  }

  bool mask(int x, int y, int z) {
    int index = x + y * mRes + z * mRes * mRes;
    return mchildMask[index];
  }

  virtual BLOCK_T *allocate(int x, int y, int z) {
    int index = x + y * mxRes + z * mxRes * myRes;
    assert(index < msize);
    if (!mchildMask[index]) {
      mchildMask[index] = 1;
      BLOCK_T *nblock = new BLOCK_T();

      //        std::cout << "  allocate: " << index << " " << msize
      //         << " - " << x << " " << y << " " << z <<  " " <<  nblock <<
      //         std::endl;

      mblocks[index] = nblock;
      return nblock;
    } else
      return mblocks[index];
  }

  virtual BLOCK_T *get(int x, int y, int z) {
    int index = x + y * mxRes + z * mxRes * myRes;

    if (!mchildMask[index]) {
      return this->allocate(x, y, z);
    };
    return mblocks[index];
  }

  virtual bool get(int x, int y, int z, BLOCK_T *block,
                   std::function<bool(coordinate_type, T, void *)> cb,
                   void *data) {
    int index = x + y * mxRes + z * mxRes * myRes;
    // std::cout << "  get: " << x << " " << y << " " << z << std::endl;
    // std::cout << "  res: " << mxRes << " " << myRes << " " << mzRes <<
    // std::endl;

    if (!mchildMask[index] && cb(center(x, y, z), mdx, data)) {
      block = this->allocate(x, y, z);
      return true;
    } else
      return false;
  }

  virtual void nearestBin(coordinate_type pos, int *bini) {
    coordinate_type binf = pos - (mcenter - 0.5 * mlengths +
                                  0.5 * coordinate_type(mdx, mdx, mdx, 0.0));
    bini[0] = binf[0] / mdx + 0.5;
    bini[1] = binf[1] / mdx + 0.5;
    bini[2] = binf[2] / mdx + 0.5;
  }

  vector<int> intersectingBins(box_type &bbox) {}

  // binning variables
  BLOCK_T *mblocks[N * N * N]; // these can be null;
  std::bitset<N * N * N> mchildMask;
  int msize;
  int mRes;
  int mxRes;
  int myRes;
  int mzRes;
  int mnumBlocks;
  T mdx;

  coordinate_type mcenter;
  coordinate_type mlengths;
};

// T                             128
// L0 = 4           32       32       32       32
// L1 = 4        8 8 8 8  8 8 8 8  8 8 8 8  8 8 8 8
// L2 = 8 11111111
// 43, 86, 103
//`35: 1, 1, 3

void testBins() {
  int logNx0 = 2;
  int logNx1 = 2;
  int logNx2 = 3;
  int N0 = pow(2, logNx0);
  int N1 = pow(2, logNx1);
  int N2 = pow(2, logNx2);

  std::cout << "---test bins---" << std::endl;
  std::cout << N0 << " " << N1 << " " << N2 << std::endl;
  std::cout << N0 << " " << N1 * N0 << " " << N2 * N1 * N0 << std::endl;

  int xi = 1753;
  int i0 = xi >> logNx2 >> logNx1;
  int r1 = xi - (i0 << logNx2 << logNx1);
  int i1 = r1 >> logNx2;
  int i2 = (r1 - (i1 << logNx2));

  std::cout << (1 << logNx0) << " " << (1 << logNx0 << logNx1) << " "
            << (1 << logNx0 << logNx1 << logNx2) << std::endl;

  std::cout << xi << std::endl;
  std::cout << i0 << " " << i1 << " " << i2 << std::endl;

  std::cout << (i2 + (i1 << logNx2) + (i0 << logNx1 << logNx2)) << " "
            << std::endl;

  std::cout << "---------------" << std::endl;
}

template <typename SPACE>
Block<SPACE, 24, std::vector<asawa::vertex<SPACE> *>> *
makeBlock(asawa::surf<SPACE> *mesh) {
  M2_TYPEDEFS;
  typedef Block<SPACE, 24, std::vector<vertex_ptr>> vertex_block;

  testBins();
  bounding_box<SPACE> bb = mesh->calc_bbox();
  coordinate_type cen = 0.5 * bb.min + 0.5 * bb.max;
  cen[3] = 1.0;
  coordinate_type dl = bb.max - bb.min;
  T mdl = max(dl[2], max(dl[0], dl[1]));

  vertex_block *block = new vertex_block(mdl, cen);
  vector<vertex_ptr> verts = mesh->get_vertices();
  std::cout << "verts! " << verts.size() << std::endl;

  for (int i = 0; i < verts.size(); ++i) {
    int b[3];
    block->nearestBin(verts[i]->coordinate(), b);
    vector<vertex_ptr> &bverts = block->get(b[0], b[1], b[2]);
    bverts.push_back(verts[i]);
  }
  return block;
};

template <typename SPACE, int L0, int L1, int L2, typename BLOCK_T>
class BlockTree
    : public Block<
          SPACE, L0,
          Block<SPACE, L1, Block<SPACE, L2, LeafBlock<SPACE, BLOCK_T>>>> {

  M2_TYPEDEFS;

public:
  typedef Block<SPACE, L0,
                Block<SPACE, L1, Block<SPACE, L2, LeafBlock<SPACE, BLOCK_T>>>>
      Root;
  typedef Block<SPACE, L1, Block<SPACE, L2, LeafBlock<SPACE, BLOCK_T>>> Block1;
  typedef Block<SPACE, L2, LeafBlock<SPACE, BLOCK_T>> Block2;
  typedef LeafBlock<SPACE, BLOCK_T> Leaf;
  typedef Eigen::Vector2i Vec2i;
  typedef Eigen::Vector3i Vec3i;

  BlockTree(T length, coordinate_type cen)
      : Root(length, cen), mLN0(log2(L0)), mLN1(log2(L1)), mLN2(log2(L2)) {

    msSize = 1 << mLN0 << mLN1 << mLN2;

    msdx = length / msSize;
  }

  Vec3i getLevels(const int &xi) {

    int i0 = xi >> mLN2 >> mLN1;
    int r1 = xi - (i0 << mLN2 << mLN1);
    int i1 = r1 >> mLN2;
    int i2 = (r1 - (i1 << mLN2));
    return Vec3i(i0, i1, i2);
  }

  virtual Leaf *topGet(int x, int y, int z) {
    // int nTotal =

    Vec3i xi = getLevels(x);
    Vec3i yi = getLevels(y);
    Vec3i zi = getLevels(z);

    // int index0 = xi[0] + yi[0]*mLN0 + zi[0]*mLN0*mLN0;
    // int index1 = xi[1] + yi[1]*mLN1 + zi[1]*mLN1*mLN1;
    // int index2 = xi[2] + yi[2]*mLN2 + zi[2]*mLN2*mLN2;
    Root &root = *this;
    Block1 *b1 = Root::get(xi[0], yi[0], zi[0]);
    b1->setCenter(Root::center(xi[0], yi[0], zi[0]));
    b1->setLength(Root::dx());
    Block2 *b2 = b1->get(xi[1], yi[1], zi[1]);
    b2->setCenter(b1->center(xi[1], yi[1], zi[1]));
    b2->setLength(b1->dx());
    Leaf *leaf = b2->get(xi[2], yi[2], zi[2]);
    return leaf;
  }

  virtual bool topGet(int x, int y, int z, Leaf *&leaf,
                      std::function<bool(coordinate_type, T, void *)> cb,
                      void *data) {
    // int nTotal =
    if (x < 0 || y < 0 || z < 0 || x >= msSize || y >= msSize || z >= msSize)
      return false;
    Vec3i xi = getLevels(x);
    Vec3i yi = getLevels(y);
    Vec3i zi = getLevels(z);

    // std::cout << xi.transpose() << std::endl;
    // std::cout << yi.transpose() << std::endl;
    // std::cout << zi.transpose() << std::endl;

    Root &root = *this;
    if (!cb(Root::center(xi[0], yi[0], zi[0]), Root::dx(), data))
      return false;
    // std::cout << " b1 " << std::endl;
    Block1 *b1 = Root::get(xi[0], yi[0], zi[0]);
    b1->setCenter(Root::center(xi[0], yi[0], zi[0]));
    b1->setLength(Root::dx());

    if (!cb(b1->center(xi[1], yi[1], zi[1]), b1->dx(), data))
      return false;
    // std::cout << " b2 " << std::endl;
    Block2 *b2 = b1->get(xi[1], yi[1], zi[1]);
    b2->setCenter(b1->center(xi[1], yi[1], zi[1]));
    b2->setLength(b1->dx());

    if (!cb(b2->center(xi[2], yi[2], zi[2]), b2->dx(), data))
      return false;
    leaf = b2->get(xi[2], yi[2], zi[2]);
    // std::cout << "dx: " << b2->dx() << std::endl;
    // std::cout << " leaf: "  << xi[2] << " " <<  yi[2] << " " <<  zi[2] << " "
    // << leaf << std::endl;
    return true;
  }

  coordinate_type decFloor(coordinate_type c, T dx) {
    // std::cout << " decFllor"<< c[0] << " " << dx << " "<< c[0]*dx <<  " " <<
    // floor(c[0]/dx) << " " << floor(c[0]/dx)*dx << std::endl;
    return coordinate_type(floor(c[0] / dx) * dx, floor(c[1] / dx) * dx,
                           floor(c[2] / dx) * dx, 1);
  }

  coordinate_type decCeil(coordinate_type c, T dx) {
    // std::cout << " decCeil"<< c[0] << " " << dx << " " << ceil(c[0]/dx) << "
    // " << ceil(c[0]/dx)*dx << std::endl;
    return coordinate_type(ceil(c[0] / dx) * dx, ceil(c[1] / dx) * dx,
                           ceil(c[2] / dx) * dx, 1);
  }

  void boundingIndices(bounding_box<SPACE> bb, int *binb, int *bine) {
    // coordinate_type min = decFloor(bb.min, msdx);
    // coordinate_type max = decCeil(bb.max, msdx);

    coordinate_type binMin =
        bb.min - (this->mcenter - 0.5 * this->mlengths -
                  0.5 * coordinate_type(msdx, msdx, msdx, 0.0));
    coordinate_type binMax =
        bb.max - (this->mcenter - 0.5 * this->mlengths -
                  0.5 * coordinate_type(msdx, msdx, msdx, 0.0));

    binb[0] = binMin[0] / msdx - 0.5;
    binb[1] = binMin[1] / msdx - 0.5;
    binb[2] = binMin[2] / msdx - 0.5;
    bine[0] = binMax[0] / msdx + 0.5;
    bine[1] = binMax[1] / msdx + 0.5;
    bine[2] = binMax[2] / msdx + 0.5;
#if 0
      std::cout << "== expand bbox == " << msdx << std::endl;
      std::cout << 1.0/msdx*binMin.transpose() << std::endl;
      std::cout << 1.0/msdx*binMax.transpose() << std::endl;
      std::cout << binb[0] << " " << binb[1] << " " << binb[2] << std::endl;
      std::cout << bine[0] << " " << bine[1] << " " << bine[2] << std::endl;
#endif
  }

  void nearestBin(coordinate_type pos, int *bini) {
    coordinate_type binf = pos - (this->mcenter - 0.5 * this->mlengths +
                                  0.5 * coordinate_type(msdx, msdx, msdx, 0.0));
    bini[0] = binf[0] / msdx + 0.5;
    bini[1] = binf[1] / msdx + 0.5;
    bini[2] = binf[2] / msdx + 0.5;
  }

  T msdx;
  int msSize;
  int mLN0, mLN1, mLN2;
};

template <typename SPACE, int L0, int L1, int L2>
BlockTree<SPACE, L0, L1, L2, std::vector<asawa::vertex<SPACE> *>> *
makeTree(asawa::surf<SPACE> *mesh) {
  M2_TYPEDEFS;
  typedef BlockTree<SPACE, L0, L1, L2, std::vector<asawa::vertex<SPACE> *>>
      vertex_tree;
  typedef typename vertex_tree::Root Root;
  typedef typename vertex_tree::Block1 Block1;
  typedef typename vertex_tree::Block2 Block2;

  bounding_box<SPACE> bb = mesh->calc_bbox();
  coordinate_type cen = 0.5 * bb.min + 0.5 * bb.max;
  cen[3] = 1.0;
  coordinate_type dl = bb.max - bb.min;
  T mdl = max(dl[2], max(dl[0], dl[1]));

  vertex_tree *block = new vertex_tree(mdl, cen);
  vector<vertex_ptr> verts = mesh->get_vertices();

  for (int i = 0; i < verts.size(); ++i) {
    int b[3];
    block->nearestBin(verts[i]->coordinate(), b);
    vector<vertex_ptr> &vbin = block->topGet(b[0], b[1], b[2]);
    vbin.push_back(verts[i]);
  }
  return block;
};

template <typename SPACE, int L0, int L1, int L2>
BlockTree<SPACE, L0, L1, L2, std::vector<asawa::face<SPACE> *>> *
// makeTriTree(asawa::surf<SPACE> *  mesh, ImmediateLines * debug,
// ImmediateLines * debug0){
makeTriTree(asawa::surf<SPACE> *mesh) {

  M2_TYPEDEFS;
  typedef BlockTree<SPACE, L0, L1, L2, std::vector<asawa::face<SPACE> *>>
      vertex_tree;
  typedef typename vertex_tree::Root Root;
  typedef typename vertex_tree::Block1 Block1;
  typedef typename vertex_tree::Block2 Block2;
  typedef typename vertex_tree::Leaf Leaf;

  bounding_box<SPACE> bb = mesh->calc_bbox();
  coordinate_type cen = 0.5 * bb.min + 0.5 * bb.max;
  cen[3] = 1.0;
  coordinate_type dl = bb.max - bb.min;
  T mdl = max(dl[2], max(dl[0], dl[1]));

  vertex_tree *block = new vertex_tree(mdl, cen);
  vector<face_ptr> faces = mesh->get_faces();

  auto tribox = [&](coordinate_type cen, T half, void *data) {
    coordinate_type *tri = static_cast<coordinate_type *>(data);
    va::tri_box<T, coordinate_type> tribox;

    coordinate_type halfc = 0.5 * coordinate_type(half, half, half, 0.0);

#if 0
      double cenf[3] = {cen[0], cen[1], cen[1]};
      double halff[3] = {half, half, half};
      double trif[3][3] = {
	{tri[0][0],tri[0][1],tri[0][2]},
	{tri[1][0],tri[1][1],tri[1][2]},
	{tri[2][0],tri[2][1],tri[2][2]}
      };

      
      bool inBox0 = (tri[0][0] > cen[0] - half && tri[0][0] < cen[0] + half &&
		     tri[0][1] > cen[1] - half && tri[0][1] < cen[1] + half &&
		     tri[0][2] > cen[2] - half && tri[0][2] < cen[2] + half);
      bool inBox1 = (tri[1][0] > cen[0] - half && tri[1][0] < cen[0] + half &&
		     tri[1][1] > cen[1] - half && tri[1][1] < cen[1] + half &&
		     tri[1][2] > cen[2] - half && tri[1][2] < cen[2] + half);
      bool inBox2 = (tri[2][0] > cen[0] - half && tri[2][0] < cen[0] + half &&
		     tri[2][1] > cen[1] - half && tri[2][1] < cen[1] + half &&
		     tri[2][2] > cen[2] - half && tri[2][2] < cen[2] + half); 
     
      std::cout << "==tribox==" << std::endl;
      std::cout << tribox.triBoxOverlap(cen,halfc,tri) << std::endl;
      std::cout << triBoxOverlap(cenf, halff, trif) << std::endl;
      std::cout << inBox0 << " " << inBox1 << " " << inBox2 << std::endl;
      std::cout << cen.transpose() << " - " << halfc.transpose() << std::endl;
      std::cout << " 0: " << tri[0].transpose() << std::endl;
      std::cout << " 1: " << tri[1].transpose() << std::endl;
      std::cout << " 2: " << tri[2].transpose() << std::endl;

#endif
    return tribox.triBoxOverlap(cen, halfc, tri);
  };

  // for(int i = 5; i < 6;  ++i ){

  std::cout << " num faces: " << faces.size() << std::endl;
  // for(int k = 0; k < 10;  ++k ){
  // int i = rand()%faces.size();
  for (int i = 0; i < faces.size(); ++i) {

    // insertDebugTri(*mesh, i, debug0);

    int binb[3];
    int bine[3];
    bounding_box<SPACE> bb = faces[i]->calc_bbox();
    block->boundingIndices(bb, binb, bine);

    /*
    std::cout << i
              << " bb - min: " << bb.min.transpose()
              << " max: "      << bb.max.transpose() << std::endl;
    std::cout << " " << binb[0] << " " << binb[1] << " " << binb[2] <<
    std::endl; std::cout << " " << bine[0] << " " << bine[1] << " " << bine[2]
    << std::endl;
    */

    coordinate_type *tri = new coordinate_type[3];
    tri[0] = faces[i]->fbegin()->coordinate();
    tri[1] = faces[i]->fbegin()->next()->coordinate();
    tri[2] = faces[i]->fbegin()->prev()->coordinate();

    for (int z = binb[2]; z < bine[2]; ++z)
      for (int y = binb[1]; y < bine[1]; ++y)
        for (int x = binb[0]; x < bine[0]; ++x) {
          Leaf *leaf = NULL;
          int b[3] = {x, y, z};
          if (block->topGet(b[0], b[1], b[2], leaf, tribox,
                            static_cast<void *>(tri))) {
            leaf->data.push_back(faces[i]);
          }
        }
    // delete tri;
  }
  return block;
};

} // namespace asawa
#endif
