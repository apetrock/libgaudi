//
//  m2Dense.h
//  Manifold
//
//  Created by John Delaney on 5/6/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

class m2Dense {
    m2Dense(){}
    ~m2Dense(){}
public:
    m2Control*          mControlMesh;
    vector<int>         mSubHash;
    vector<Vec<3,T> >   mVertexPool;
};