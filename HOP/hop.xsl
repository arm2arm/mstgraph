<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/TR/WD-xsl">
  <xsl:template match="/">
    <html>
      <body>
        <table border="2" bgcolor="yellow">
          <tr>
            <th>TEST</th>
            <th>VALUE</th>
          </tr>
          <xsl:for-each select="OPTIONS/test">
            <tr>
              <td>
                <xsl:value-of select="path"/>
              </td>
              <td>
                <xsl:value-of select="file"/>
              </td>
            </tr>
          </xsl:for-each>
        </table>
      </body>
    </html>
  </xsl:template>
</xsl:stylesheet>
