using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Text.RegularExpressions;
using System.Linq;


namespace MinSunVolume
{
    public class MinSunVolumeComponent : GH_Component
    {
        public MinSunVolumeComponent()
          : base("MinSunVolume", "SunVol",
              "Description",
              "Category", "Subcategory")
        {
        }
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("地理信息库", "地理信息库", "", GH_ParamAccess.list);    //0
            pManager.AddTextParameter("地点", "地点", "", GH_ParamAccess.item);            //1
            pManager.AddIntervalParameter("计算时间段", "时间段", "", GH_ParamAccess.item);         //2
            pManager.AddNumberParameter("计算精度", "精度", "", GH_ParamAccess.item);          //3
            pManager.AddTextParameter("日期", "日期", "", GH_ParamAccess.item);                   //4
            pManager.AddPointParameter("被遮挡点", "被遮挡点", "", GH_ParamAccess.list);           //5
            pManager[5].Optional = true;
            pManager.AddBrepParameter("被遮挡面", "被遮挡面", "", GH_ParamAccess.list);          //6
            pManager[6].Optional = true;
            pManager.AddNumberParameter("被遮挡面窗台计算高度", "窗台高度", "", GH_ParamAccess.item, 900); //7
            pManager[7].Optional = true;
            pManager.AddNumberParameter("被遮挡面采样间距", "采样间距", "", GH_ParamAccess.item, 5000);    //8
            pManager[8].Optional = true;
            pManager.AddMeshParameter("其他遮挡物", "其他遮挡物", "", GH_ParamAccess.list);     //9
            pManager.AddCurveParameter("形体生成范围曲线", "范围曲线", "", GH_ParamAccess.item);         //10
            pManager.AddNumberParameter("生成单元格大小", "单元格大小", "", GH_ParamAccess.item);         //11
            pManager.AddNumberParameter("最大建筑高度", "最大建筑高度", "", GH_ParamAccess.item); //12
        }
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddBoxParameter("最大无遮挡建设体量", "无遮挡建设体量", "", GH_ParamAccess.list);
            //pManager.AddVectorParameter("_Vectors", "_Vec", "", GH_ParamAccess.list);
            //pManager.AddPointParameter("_Grids", "_Gri", "", GH_ParamAccess.list);
            pManager.AddPointParameter("窗位计算点", "窗位点", "", GH_ParamAccess.list);
        }
        Curve _SiteCrv;
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region 参数初始化与数据读取
            List<Point3d> _TargetPts = new List<Point3d>();
            List<Brep> _TargetFaces = new List<Brep>();
            double _WindowElevation = 0;
            double _FacePtDistance = 1000;
            List<string> _LocationDataBase = new List<string>();
            string _Location = "";
            string _Day = "";
            Interval _TimeSpan = new Interval();
            double _StepAccu = 0;
            double _Size = 0;
            List<Mesh> _OcclusionVolumes = new List<Mesh>();
            double _MaxBuildingHeight = 200000;

            DA.GetData("形体生成范围曲线", ref _SiteCrv);
            DA.GetDataList("被遮挡点", _TargetPts);
            DA.GetDataList("被遮挡面", _TargetFaces);
            DA.GetData("被遮挡面窗台计算高度", ref _WindowElevation);
            DA.GetData("被遮挡面采样间距", ref _FacePtDistance);
            DA.GetDataList("地理信息库", _LocationDataBase);
            DA.GetData("地点", ref _Location);
            DA.GetData("日期", ref _Day);
            DA.GetData("计算时间段", ref _TimeSpan);
            DA.GetData("计算精度", ref _StepAccu);
            DA.GetData("生成单元格大小", ref _Size);
            DA.GetDataList("其他遮挡物", _OcclusionVolumes);
            DA.GetData("最大建筑高度", ref _MaxBuildingHeight);
            #endregion


            #region 启动检查
            //如果场地曲线非法，或曲线未封闭则退出
            if (_SiteCrv == null || !_SiteCrv.IsValid || !_SiteCrv.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Site is Invalid or not Closed");   //*
                return;
            }
            //如果曲线不是平面曲线,退出
            Plane tempPl;
            if (!_SiteCrv.TryGetPlane(out tempPl))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Site isn't a Planar Curve");   //*
                return;
            }
            #endregion

            #region 1\2.获取计算点
            //可计算的点包含两部分:
            //1.直接给定的计算点,该部分计算位置不补偿窗台高度
            //2.输入给定的计算面,通过对计算面在指定高度及指定间距采样得到的计算点
            List<Point3d> calPts = new List<Point3d>();  //建立计算点容器
            List<double> calMaxHeight = new List<double>();
            double tol = DocumentAngleTolerance();  //*

            //1.当指定计算点存在时,直接添加至计算点容器
            //2.当给定的计算面存在时,进行计算点采样
            if (_TargetFaces.Count != 0)
            {
                //2.1对计算面的合法性进行评估,舍弃非法的计算面
                List<Brep> targetFaces = new List<Brep>();
                foreach (Brep b in _TargetFaces)
                {
                    if (b.IsValid) targetFaces.Add(b);
                }

                if (targetFaces.Count > 0)
                {
                    //2.2获取底边标高
                    foreach (Brep b in targetFaces)
                    {
                        Curve[] cutCrvs;  //与指定高度切割后的曲线
                        Point3d[] tempPts; //函数的临时容器，无作用

                        BoundingBox box = b.GetBoundingBox(true);
                        double zMin = box.Min.Z;  //曲面底标高
                        zMin += _WindowElevation; //依照窗台高度偏移后标高
                        double zMax = box.Max.Z - zMin;

                        //2.3建立高度切割平面
                        Plane pl = Plane.WorldXY;
                        pl.Translate(new Vector3d(0, 0, zMin));

                        //2.4生成高度切割线
                        Rhino.Geometry.Intersect.Intersection.BrepPlane(
                          b, //待分割的计算面
                          pl, //高度切割平面
                          tol, //模型公差
                          out cutCrvs, //切割曲线结果容器
                          out tempPts); //切割点结果临时容器

                        //2.5对切割曲线依照指定间距划分
                        try
                        {
                            foreach (Curve c in cutCrvs)
                            {
                                double[] t = c.DivideByLength(_FacePtDistance, true); //依照距离生成计算点的T值
                                foreach (double d in t)
                                {
                                    calPts.Add(c.PointAt(d));  //根据T值得到具体计算点
                                    calMaxHeight.Add(zMax);
                                }
                            }
                        }
                        catch
                        { }
                    }
                }
            }
            #endregion

            #region 3.获取地理信息并得出计算光线
            //3.1解析地理位置信息
            double attitude, longttitude;  //计算地点的纬度\精度参数
            string loc;                   //计算点的地理名称
            GetLocationInfo(_Location, _LocationDataBase, out attitude, out longttitude, out loc);  //地理信息获取方法

            double dayLatti;  //计算地点计算时间的赤纬角参数
            GetDayLattitude(_Day, out dayLatti);  //加载赤纬计算方法

            int timeCount;  //计算光线总数参数
            Vector3d[] rays = GenearteRay(
              attitude, //纬度参数
              dayLatti, //赤纬参数
              _TimeSpan.Min, _TimeSpan.Max, //计算时间区间
              _StepAccu, //计算精度
              out timeCount//计算光线总数
              );
            //Debug 计算光线线输出
            //_Vectors = rays;  //*
            #endregion

            #region 4.场地曲线划分
            int u;  //横向分割数量
            int v;  //纵向风格数量

            Point3d[] Grids = GenerateGrids(//加载场地曲线分割方法
              _SiteCrv, //场地曲线
              _Size, //划分单元尺寸
              out u, //横向分割数量
              out v//纵向风格数量
              );
            #endregion

            #region 周边遮挡预处理
            List<Mesh> OcclusionVolumes = new List<Mesh>();
            foreach (Mesh m in _OcclusionVolumes)
            {
                if (m == null) continue;
                if (m.IsValid) OcclusionVolumes.Add(m);
            }
            #endregion

            #region 5.计算每个采样点的不遮挡最小高度
            //结果容器
            List<double> timings = new List<double>();
            List<double> maxHeights = new List<double>();
            List<Line> debugs = new List<Line>();

            //5.1设置最大建筑高度
            double heightLimit = 200000;
            if (_MaxBuildingHeight > double.Epsilon) heightLimit = _MaxBuildingHeight;

            for (int i = 0; i < Grids.Length; i++)
            {
                //5.2核心解算
                double heightG = heightLimit;
                double time;        //计算时点
                double height;      //解算高度
                for (int j = 0; j < calPts.Count; j++)
                {
                    Vector3d CordiVector = calPts[j] - Grids[i];  //遮挡点与计算点的指向向量
                    Vector3d ptSunRay;  //采样计算光线


                    if (Grids[i] != Point3d.Unset) //检查计算点的合法性
                    {
                        ptSunRay = SunRayByCordinateVector(rays, CordiVector, out time); //解算采样计算光线
                        Vector3d occRay = ptSunRay;
                        double occHeightOffset = 0;
                        occRay.Reverse();

                        bool Occlued = false;
                        if (OcclusionVolumes.Count != 0)
                        {
                            Occlued = OcclusionTest(occRay, calPts[j], Grids[i], OcclusionVolumes, tol, out occHeightOffset);
                        }

                        height = MaxHeight(calPts[j], Grids[i], ptSunRay);//解算计算点最大高度*

                        if (Occlued && occHeightOffset <= calMaxHeight[j])
                        {
                            height += occHeightOffset;
                        }
                        else if (Occlued && occHeightOffset > calMaxHeight[j])
                        {
                            height = heightLimit;
                        }
                    }
                    else  //非法计算点的处理
                    {
                        ptSunRay = Vector3d.Unset;
                        height = double.PositiveInfinity;
                        time = double.PositiveInfinity;
                    }

                    timings.Add(time);  //记录解算时间
                    if (height < heightG) heightG = height;  //记录最小高度
                }
                maxHeights.Add(heightG);  //记录最小高度结果
            }

            List<Point3d> LGrids = Grids.ToList();  //Grids数组转换为列表
            List<BoundingBox> results = GenerateData(LGrids, maxHeights, u, v); //调用单元体生成方法
            #endregion

            #region 6.结果输出
            DA.SetDataList("最大无遮挡建设体量", results);
            //DA.SetDataList("_Grids", Grids);
            DA.SetDataList("窗位计算点", calPts);
            //DA.SetDataList("_Vectors", rays);

            #endregion
        }
        #region 获取地理位置方法，返回是否成功
        public static bool GetLocationInfo(
          string regStr, //地址检索条件（正则语句）
          List<string> dataBase, //地理信息数据
          out double attitude, //out 纬度
          out double longttitude, //out 经度
          out string location// 检索到的地点
          )
        {
            string dataStr = "";
            bool findSucecess = false;
            for (int i = 0; i < dataBase.Count; i++)
            {
                string regexStr = regStr;   //设置正则表达式语句
                Regex reg = new Regex(regexStr, RegexOptions.None); //创建正则表达式类
                Match matchClass = reg.Match(dataBase[i]);

                if (matchClass.Success)
                {
                    findSucecess = true;
                    dataStr = dataBase[i];
                    break;
                }
            }
            string[] part = dataStr.Split('|');
            location = part[0] + ":" + part[1];
            longttitude = Convert.ToDouble(part[2]);
            attitude = Convert.ToDouble(part[3]);
            return findSucecess;
        }
        #endregion

        #region 获取计算日期的赤纬角方法，返回是否成功
        public static bool GetDayLattitude(
          string day, //计算日期，考虑到国标的实际使用，此处仅提供“大寒日”与“冬至日”作为输入口（未来改为Enum）
          out double dayLattitude//out 计算日期的赤纬角
          )
        {
            bool sucecess = false;
            Dictionary<string, double> dayValue = new Dictionary<string, double> { { "大寒日", -20.283 }, { "冬至日", -23.45 } };
            sucecess = dayValue.TryGetValue(day, out dayLattitude);
            return sucecess;
        }
        #endregion

        #region 生成计算光线方法，返回光线向量数组
        public static Vector3d[] GenearteRay(
          double lattitude, //纬度参数
          double dayLatti, //赤纬参数
          double stTime, //计算起始时间（小时)
          double endTime, //计算结束时间（小时）
          double accur, //计算精度（分钟）
          out int timeCount//out 生成光线总数
          )
        {
            double attiRad = lattitude / 180 * Math.PI;
            double dayRad = dayLatti / 180 * Math.PI;

            int oTimeCount = (int)Math.Ceiling((endTime - stTime) * 60 / accur);

            Vector3d[] rays = new Vector3d[oTimeCount + 1];
            double[] times = new double[oTimeCount + 1];

            for (int t = 0; t <= oTimeCount; t++)
            {
                double time1 = t * accur + stTime * 60;
                double angTime = (15 * (time1 / 60 - 12)) / 180 * Math.PI;

                double sinH = Math.Sin(attiRad) * Math.Sin(dayRad) + Math.Cos(attiRad) * Math.Cos(dayRad) * Math.Cos(angTime);
                double angH = Math.Asin(sinH);

                double cosA = (sinH * Math.Sin(attiRad) - Math.Sin(dayRad)) / (Math.Cos(angH) * Math.Cos(attiRad));
                double angA = Math.Acos(cosA);
                if (angTime > double.Epsilon)
                {
                    angA = -angA;
                }

                Vector3d ray = Vector3d.YAxis;
                ray.Rotate(-angH, Vector3d.XAxis);
                ray.Rotate(angA, Vector3d.ZAxis);
                ray.Unitize();

                rays[t] = ray;
            }
            timeCount = oTimeCount;
            return rays;
        }
        #endregion

        #region 根据计算点与被遮挡点关系，确定其计算日时刻及计算光线（插值计算）
        public static Vector3d SunRayByCordinateVector(
          Vector3d[] rays, //基准光线（向量）序列数组
          Vector3d testRay, //平面指向向量(被遮挡点指向计算点)
          out double timing//out 计算时刻
          )
        {
            timing = double.NaN;
            double[] angleRaysToX = new double[rays.Length];

            for (int i = 0; i < rays.Length; i++)
            {
                Vector3d tempRay = new Vector3d(rays[i].X, rays[i].Y, 0);
                double angleTemp = Vector3d.VectorAngle(tempRay, Vector3d.XAxis, Plane.WorldXY);
                angleRaysToX[i] = angleTemp;
            }

            Vector3d planarTestRay = new Vector3d(testRay.X, testRay.Y, 0);
            double angleTest = Vector3d.VectorAngle(planarTestRay, Vector3d.XAxis, Plane.WorldXY);

            if (angleTest < angleRaysToX[0] || angleTest > angleRaysToX[angleRaysToX.Length - 1])
            {
                return Vector3d.Unset;
            }

            Vector3d VecTemp = Vector3d.Unset;
            for (int i = 1; i < angleRaysToX.Length; i++)
            {
                if (angleTest == angleRaysToX[i] || (angleTest > angleRaysToX[i - 1] && angleTest < angleRaysToX[i]))
                {
                    double t = (angleTest - angleRaysToX[i - 1]) / (angleRaysToX[i] - angleRaysToX[i - 1]);
                    VecTemp = (rays[i] - rays[i - 1]) * t + rays[i - 1];
                    timing = t + i - 1;
                    break;
                }
            }
            return VecTemp;
        }
        #endregion

        #region 遮挡物遮挡测试，返回是否被遮挡
        public static bool OcclusionTest(Vector3d rayVec, Point3d orient, Point3d targetPt, List<Mesh> Volumes, double tolerence, out double occHeight)
        {
            bool Occlued = false;
            rayVec.Unitize();
            Vector3d flatRay = new Vector3d(rayVec.X, rayVec.Y, 0);
            double ratio = flatRay.Length / rayVec.Length;
            Vector3d cordiRay = targetPt - orient;
            cordiRay.Z = 0;
            double ptDis = cordiRay.Length;
            orient = orient + (tolerence * rayVec);
            double tempHeightOffset = 0;

            foreach (Mesh m in Volumes)
            {
                Ray3d testRay0 = new Ray3d(orient, rayVec);
                double dis0 = Rhino.Geometry.Intersect.Intersection.MeshRay(m, testRay0);
                if (dis0 >= 0 && dis0 * ratio <= ptDis)
                {
                    Occlued = true;

                    Point3d tempPt = orient;
                    bool tempOcc = true;
                    double offsetVal = 0;
                    double offsetStep = 1000;
                    while (tempOcc)
                    {
                        Ray3d testRay = new Ray3d(tempPt, rayVec);
                        double dis = Rhino.Geometry.Intersect.Intersection.MeshRay(m, testRay);
                        if (dis >= 0 && dis * ratio <= ptDis)
                        {
                            tempPt = tempPt + new Vector3d(0, 0, offsetStep);
                            offsetVal += offsetStep;
                        }
                        else
                        {
                            tempOcc = false;
                        }
                        if (offsetVal > tempHeightOffset)
                        {
                            tempHeightOffset = offsetVal;
                        }
                    }
                }
            }
            occHeight = tempHeightOffset;
            return Occlued;
        }
        #endregion

        #region 根据被遮挡点、计算点、阳光向量解算最大建设高度
        public static double MaxHeight(
          Point3d orient, //被遮挡点
          Point3d testPoint, //计算点
          Vector3d sunRay//计算太阳向量
          )
        {
            Vector3d dir = orient - testPoint;
            double zDif = dir.Z;
            Vector3d planarSunRay = new Vector3d(sunRay.X, sunRay.Y, 0);
            double t = sunRay.Z / planarSunRay.Length;
            double height = t * (new Vector3d(dir.X, dir.Y, 0)).Length;
            double tempHeight = zDif - height;
            return tempHeight;
        }
        #endregion

        #region 基地曲线拆分成网格阵列
        public static Point3d[] GenerateGrids(
          Curve boundaryCrv, //基地曲线
          double size, //网格单元尺寸(mm)
          out int uDiv, //out 横向划分单元数
          out int vDiv//out 纵向划分单元数
          )
        {
            Plane pl;
            uDiv = 0;
            vDiv = 0;
            if (!boundaryCrv.TryGetPlane(out pl))
            {
                return new Point3d[] { Point3d.Unset };
            }
            BoundingBox box = boundaryCrv.GetBoundingBox(true);
            List<Point3d> pts = new List<Point3d>();
            double uDis = box.Max.X - box.Min.X;
            double vDis = box.Max.Y - box.Min.Y;

            uDiv = Convert.ToInt32(Math.Ceiling(uDis / size));
            vDiv = Convert.ToInt32(Math.Ceiling(vDis / size));
            for (int i = 0; i < vDiv; i++)
            {
                for (int j = 0; j < uDiv; j++)
                {
                    Point3d tempPt = new Point3d(box.Min.X + j * size, box.Min.Y + i * size, 0);
                    tempPt.Transform(Transform.ProjectAlong(pl, Vector3d.ZAxis));
                    var a = boundaryCrv.Contains(tempPt, pl, 0.1);
                    if (a == PointContainment.Inside || a == PointContainment.Coincident)
                    {
                        pts.Add(tempPt);
                    }
                    else
                    {
                        pts.Add(Point3d.Unset);
                    }
                }
            }
            return pts.ToArray();
        }
        #endregion

        #region 依据计算结果，建立单元体量
        public static List<BoundingBox> GenerateData(
          List<Point3d> pts, //点列表
          List<double> heights, //最大建设高度列表
          int u, //横向划分数量
          int v
          )
        {
            List<BoundingBox> volumes = new List<BoundingBox>();

            if (pts.Count != heights.Count || pts.Count != u * v)
            {
                return volumes;
            }

            for (int i = 0; i < pts.Count - u - 1; i++)
            {
                if (i % u == u - 1) continue;
                Point3d a = pts[i];
                Point3d b = pts[i + 1];
                Point3d c = pts[i + u];
                Point3d d = pts[i + u + 1];

                if (a == Point3d.Unset || b == Point3d.Unset || c == Point3d.Unset || d == Point3d.Unset)
                {
                    continue;
                }
                Point3d e = a + new Vector3d(0, 0, heights[i]);
                Point3d f = b + new Vector3d(0, 0, heights[i + 1]);
                Point3d g = c + new Vector3d(0, 0, heights[i + u]);
                Point3d h = d + new Vector3d(0, 0, heights[i + u + 1]);

                double bottomHeight = Math.Min(Math.Min(a.Z, b.Z), Math.Min(c.Z, d.Z));
                double topHeihgt = Math.Min(Math.Min(e.Z, f.Z), Math.Min(g.Z, h.Z));

                a.Z = bottomHeight;
                h.Z = topHeihgt;

                BoundingBox box = new BoundingBox(a, h);
                volumes.Add(box);
            }
            return volumes;
        }
        #endregion
        protected override System.Drawing.Bitmap Icon => null;

        public override Guid ComponentGuid
        {
            get { return new Guid("62d48c12-3361-44a7-9feb-7b91342f0aeb"); }
        }
    }
}
