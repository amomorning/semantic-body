using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.IO;
using System;

public class RotateTemp : MonoBehaviour
{
    public GameObject AVE;
    Vector3[] vert;
    Mesh mesh;
    Slider[] mySliders;
    float[] curVa = new float[26];

    void Start()
    {
        Debug.Log("RotateTemp Start");
        Debug.Log("?????? : " + AVE.name);
        AVE.gameObject.transform.Rotate(new Vector3(270, 315, 0));
        mesh = GetComponent<MeshFilter>().mesh;
        vert = mesh.vertices;
        float min = 0;
        for (var i = 0; i < vert.Length; ++i)
        {
            float tmp = vert[i].z;
            if (tmp < min) min = tmp;
        }


        readFile(Application.dataPath + "/Resources/measure.txt");
        Debug.Log(curVa.Length + " " + curVa[0]);
        for (int i = 0; i < curVa.Length; ++i)
        {
            string name = "Slider " + (i + 1);
           // Debug.Log(name);
            GameObject.Find(name).GetComponent<Slider>().value = curVa[i];
        }

        AVE.gameObject.transform.position += new Vector3(0.0f, -1.2f - min, 0.0f);
        Debug.Log("ixxxxx");
        Debug.Log(AVE.name);
        //AVE.gameObject.SetActive(false);
    }

    void Update()
    {
        if (Input.GetMouseButton(0))
        {
            AVE.gameObject.transform.Rotate(new Vector3(0, 0, -Input.GetAxis("Mouse X") * Time.deltaTime * 300));
        }
    }

    void readFile(string fileName)
    {
        string[] strs = File.ReadAllLines(fileName);
        for (int i = 0; i < strs.Length; i++)
        {
            curVa[i] = (float)Convert.ToDouble(strs[i]);
        }
    }
}
